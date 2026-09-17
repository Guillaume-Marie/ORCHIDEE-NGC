! =================================================================================================================================
! MODULE 	: stomate_data
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "stomate_data" module defines the values about the PFT parameters. It will print
!! the values of the parameters for STOMATE in the standard outputs. 
!!
!!\n DESCRIPTION: None 
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months))
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_data.f90 $
!! $Date: 2026-01-06 18:52:43 +0100 (mar. 06 janv. 2026) $
!! $Revision: 9290 $
!! \n
!_ ================================================================================================================================

MODULE stomate_data

  ! modules used:

  USE constantes
  USE time, ONLY : one_day, dt_sechiba, one_year
  USE pft_parameters
  USE dynamic_parameters
  USE defprec
  USE function_library, ONLY: get_printlev
  USE ioipsl_para, ONLY : ipslerr_p
  USE xios_orchidee, ONLY : xios_orchidee_send_field
  USE mod_orchidee_para, ONLY : allreduce_sum, is_root_prc, nbp_glo, numout

  IMPLICIT NONE

  INTEGER(i_std), SAVE, PRIVATE                :: printlev_loc   !! local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: hori_index     !! Move to Horizontal indices
!$OMP THREADPRIVATE(hori_index)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: horipft_index  !! Horizontal + PFT indices
!$OMP THREADPRIVATE(horipft_index)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: horican_index  !! Horizontal + canopy levels indices
!$OMP THREADPRIVATE(horican_index)
  
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: horicut_index  !! Horizontal + cut times indices
!$OMP THREADPRIVATE(horicut_index)

  !
  ! Land cover change
  !
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_s_index   !! Horizontal + P short indices
!$OMP THREADPRIVATE(horip_s_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_m_index   !! Horizontal + P medium indices
!$OMP THREADPRIVATE(horip_m_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_l_index   !! Horizontal + P long indice
!$OMP THREADPRIVATE(horip_l_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_ss_index  !! Horizontal + P short indices
!$OMP THREADPRIVATE(horip_ss_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_mm_index  !! Horizontal + P medium indices
!$OMP THREADPRIVATE(horip_mm_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip_ll_index  !! Horizontal + P long indices
!$OMP THREADPRIVATE(horip_ll_index)

  INTEGER(i_std),SAVE :: itime                 !! time step
!$OMP THREADPRIVATE(itime)
  INTEGER(i_std),SAVE :: hist_id_stomate       !! STOMATE history file ID
!$OMP THREADPRIVATE(hist_id_stomate)
  INTEGER(i_std),SAVE :: hist_id_stomate_IPCC  !! STOMATE history file ID for IPCC output
!$OMP THREADPRIVATE(hist_id_stomate_IPCC)
  INTEGER(i_std),SAVE :: rest_id_stomate       !! STOMATE restart file ID
!$OMP THREADPRIVATE(rest_id_stomate)

  REAL(r_std),PARAMETER :: adapted_crit = 1. - ( 1. / euler ) !! critical value for being adapted (1-1/e) (unitless)
  REAL(r_std),PARAMETER :: regenerate_crit = 1. / euler       !! critical value for being regenerative (1/e) (unitless)

  !
  ! Guillaume M. -- Fragmentation / Average Edge Distance (Bowring et al. 2024).
  ! Two AED variables coexist because SPITFIRE and LIGHT use different denominators:
  ! AED_fire uses area_burnable_veg (trees + grass that can burn), AED_light uses
  ! area_forest (is_tree only). Both share edge_length (m) as their common input,
  ! read annually from Edge_length.nc by spitfire_annual_input.
  !
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: AED_fire    !! Avg edge distance for SPITFIRE (m)
!$OMP THREADPRIVATE(AED_fire)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: AED_light   !! Avg edge distance for canopy light effect (m)
!$OMP THREADPRIVATE(AED_light)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: area_forest !! Forest area per grid cell (ha, is_tree only)
!$OMP THREADPRIVATE(area_forest)

  !
  ! Guillaume M. -- Dynamic edge_length (AED_FEEDBACK). edge_length_dyn evolves each
  ! timestep from the perturbation accumulators (dA_fire/storm/pest, m² per dt_stomate),
  ! relaxing toward the continuous-forest baseline aed_edge_ref_frac·edge_length_ref
  ! (default 0). See update_edge_length and MODULE_DESIGN_AED_FEEDBACK.md.
  !
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: edge_length_dyn !! dynamic edge length (m)
!$OMP THREADPRIVATE(edge_length_dyn)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: edge_length_ref !! EFDA reference edge length: validation reference & static edge when feedback off (m)
!$OMP THREADPRIVATE(edge_length_ref)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: dA_fire   !! area perturbed by fire over dt_stomate (m²)
!$OMP THREADPRIVATE(dA_fire)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: dA_storm  !! area perturbed by storm over dt_stomate (m²)
!$OMP THREADPRIVATE(dA_storm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: dA_pest   !! area perturbed by pest over dt_stomate (m²)
!$OMP THREADPRIVATE(dA_pest)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: dA_harvest !! area harvested or LCC-cut over dt_stomate (m²)
!$OMP THREADPRIVATE(dA_harvest)

  !
  ! Guillaume M. -- AED_REGROWTH class-0 conveyor: stand-replaced area held as annual
  ! cohorts per agent, ageing one bin per year and graduating at n_recruit years
  ! (canopy reclosure). Under OK_EDGE_FROM_AGE_CLASS=y it replaces the tau_rec
  ! relaxation of edge_length_dyn by a mechanistic recovery and yields frac_class0.
  ! veget_max is NOT modified. See MODULE_DESIGN_AED_REGROWTH.md.
  !
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: class0_area      !! disturbed-area cohorts (npts,nagent,n_recruit_max) (m²)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:) :: class0_veg       !! MATURITY_CONVEYOR : aire en regeneration (npts,nvmap,n_class0),
                                                                    !! Guillaume M. -- regenerating area as a grid-cell FRACTION, indexed by
                                                                    !! forest GROUP -- not by agent like class0_area -- so that the area is
                                                                    !! returned to the right species after n_class0 years. Both conveyors
                                                                    !! coexist: class0_area feeds the edge, class0_veg carries actual area.
!$OMP THREADPRIVATE(class0_veg)
!$OMP THREADPRIVATE(class0_area)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)     :: class0_recovered !! area graduating back to forest this year (m²/yr)
!$OMP THREADPRIVATE(class0_recovered)
  ! Guillaume M. -- MATURITY_CONVEYOR: grid-cell fraction returned by the class-0 conveyor
  ! to the age-class-1 slot, per slot. Published by stomate_lpj and consumed later in the
  ! SAME annual pass by land_cover_change_main, so call order replaces restart persistence.
  ! Needed because loss_gain carries only ONE scalar per slot: netting feed against
  ! restitution erases the gain event that triggers replanting. This variable restores it.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)   :: class0_restitute !! fraction restituee en classe 1 cette annee (-)
!$OMP THREADPRIVATE(class0_restitute)
  ! Guillaume M. -- AED_REGROWTH V2: edge derived from the young age-class area A_open
  ! (bounded by forest area), agent mix from a leaky disturbance-flux integrator. Replaces
  ! the unbounded class-0 conveyor. See MODULE_DESIGN_AED_REGROWTH_V2.md.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)   :: agent_flux       !! leaky-integrated dA per agent (npts,nagent) for the v2 agent mix
!$OMP THREADPRIVATE(agent_flux)
  REAL(r_std), SAVE         :: class0_year_clock_s = 0.0_r_std      !! seconds elapsed in the current conveyor year
!$OMP THREADPRIVATE(class0_year_clock_s)

  !
  ! Guillaume M. -- AED_SPINUP: spatial α lookup and spinup state.
  ! See MODULE_DESIGN_AED_SPINUP.md.
  !
  INTEGER(i_std), PARAMETER :: nagent       = 4                  !! disturbance agents (fire, storm, pest, harvest)
  INTEGER(i_std), PARAMETER :: IDX_FIRE     = 1
  INTEGER(i_std), PARAMETER :: IDX_STORM    = 2
  INTEGER(i_std), PARAMETER :: IDX_PEST     = 3
  INTEGER(i_std), PARAMETER :: IDX_HARVEST  = 4
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: mi_frac         !! fractions de classe d'intensite (npts, nmiclass), lues annuellement
!$OMP THREADPRIVATE(mi_frac)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: rotation_ref_map !! ROTATION_REF_FILE: reference rotation per point and PFT (yr),
                                                                  !! <= 0 = namelist ROTATION_REF(pft). Allocated only when the file is read.
!$OMP THREADPRIVATE(rotation_ref_map)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: road_length_map    !! FRAG_3TERMS macro: longueur de ROUTE par maille (m), statique.
                                                                  !! Guillaume M. -- NOT an edge length: update_edge_length derives
                                                                  !! EL_macro = 2 * RL * road_grip_correction from it.
!$OMP THREADPRIVATE(road_length_map)

  ! Guillaume M. -- FRAG_3TERMS permanent: lake, coast and river edge, read from a map.
  ! The meso term is built from veget_max on area_land, so its non-forest side is
  ! TERRESTRIAL ONLY: no channel exists by which water could cut forest. Static, like roads.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: edge_hydro_map     !! FRAG_3TERMS permanent: lisiere foret/eau et foret/riviere par maille (m), statique.
!$OMP THREADPRIVATE(edge_hydro_map)
  !! Guillaume M. -- AED_EDGE_RDI_WEIGHT: relative density index N/Nmax(D) per slot,
  !! PRODUCED by stomate_lpj_vegetation, where it is already computed, and CONSUMED by
  !! update_edge_length; stomate.f90 guarantees that order. rdi was recomputed locally in
  !! five modules -- this adds a shared read, not a sixth computation.
  !!
  !! Guillaume M. -- dia_factor scales the clearcut diameter with management intensity.
  !! It multiplies largest_tree_dia AND, through it, the age-class bounds that are
  !! fractions of it: the whole maturity scale of the pixel is rescaled.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: dia_factor  !! Clearcut-diameter multiplier per pixel (-), 1 = unchanged
!$OMP THREADPRIVATE(dia_factor)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: rdi_stand       !! Relative density index N/Nmax(D) per (point, PFT slot) (-)
!$OMP THREADPRIVATE(rdi_stand)
  ! Guillaume M. -- MAT_AGE_ENTRY: leaky mean of the age of the area that enters each slot on
  ! an upward class transfer (yr), -1 = never fed. A cohort attribute carried by the moved
  ! parcel, not a slot velocity. PERSISTED (libIGCM restarts the binary every period).
  ! Feeds the online a3 of the terminal setpoint and the age-based R_growth.
  ! See design/MODULE_DESIGN_MAT_AGE_ENTRY.md.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: mat_age_entry   !! Age of the area entering the slot, leaky mean (yr); <= 0 = no sample
!$OMP THREADPRIVATE(mat_age_entry)
  ! Guillaume M. -- ROTATION_GROWTH: the rotation the STAND can actually deliver, against the
  ! published itinerary. Derived from mat_age_entry once a year, persisted because it is only
  ! recomputed on 31 December. Blended with a weight w in set_management_intensity.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: rotation_growth !! Rotation derived from realised growth (yr), <=0 = undefined
!$OMP THREADPRIVATE(rotation_growth)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: clearcut_area     !! taille de parcelle de REALISATION (m2). ATTENTION : geometrie de lisiere
                                                                 !! Guillaume M. -- edge geometry ONLY: it does NOT drive what the model
                                                                 !! cuts (a cut empties the whole slot). It must never reach
                                                                 !! sapiens_forestry nor the carbon code (tested property).
!$OMP THREADPRIVATE(clearcut_area)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:) :: alpha_harvest !! coefficient de production de lisiere du harvest (m/m2)
!$OMP THREADPRIVATE(alpha_harvest)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: target_rotation_age  !! rotation cible (npts,nvm) (an) — SEUL canal agissant sur la dynamique.
                                                                          !! Guillaume M. -- PER PFT: a grid cell is a mosaic, and a single scalar
                                                                          !! cannot carry 10 yr (eucalypt) and 150 yr (oak) at the same time.
!$OMP THREADPRIVATE(target_rotation_age)
  ! Guillaume M. -- PROGRESSIVE_HARVEST: slot (ipts,ivm) area fraction actually clearcut
  ! this year by the phased clearcut (0 = no phased clearcut on that slot). Written by
  ! sapiens_forestry_main, read by sapiens_kill (harvest_area scaling) then by
  ! age_class_distr (slot split), and reset at each annual call of sapiens_forestry_main.
  ! Allocated only under ok_progressive_harvest -> test ALLOCATED before use.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: ph_cut_frac  !! fraction d'aire du slot coupee (-)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: rdi_thin_limit_diag !! RDI target actually handed to thinning, AFTER the band override (-)
                                                                    !! Guillaume M. -- RDI_TARGET_LOWER is recomputed for output without that
                                                                    !! override, so it does NOT report what management used. undef = the
                                                                    !! thinning decision point was never reached for this slot.
!$OMP THREADPRIVATE(rdi_thin_limit_diag)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: ac12_r_diag  !! MATURITY_CONVEYOR : rapport de hauteur h1/h_ref, DIAGNOSTIC (-)
!$OMP THREADPRIVATE(ac12_r_diag)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:) :: ac12_freal_diag !! MATURITY_CONVEYOR : area fraction of the donor REALLY moved
                                                                 !! this year, read AFTER the veget_max write, DIAGNOSTIC (-).
                                                                 !! Guillaume M. -- the decided rate is clipped by the guards
                                                                 !! before the write; only the realised flux is trustworthy.
!$OMP THREADPRIVATE(ac12_freal_diag)
!$OMP THREADPRIVATE(ph_cut_frac)


  ! private & public routines

  PUBLIC data, calculate_aed, update_edge_length, set_management_intensity

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: data
!!
!>\BRIEF         This routine defines the values of the PFT parameters. It will print the values of the parameters for STOMATE
!!               in the standard outputs of ORCHIDEE. 
!!
!! DESCRIPTION : This routine defines PFT parameters. It initializes the pheno_crit structure by tabulated parameters.\n
!!               Some initializations are done for parameters. The SLA is calculated according *to* Reich et al (1992).\n
!!               Another formulation by Reich et al(1997) could be used for the computation of the SLA.
!!               The geographical coordinates might be used for defining some additional parameters
!!               (e.g. frequency of anthropogenic fires, irrigation of agricultural surfaces, etc.). \n
!!               For the moment, this possibility is not used. \n
!!               The specifc leaf area (SLA) is calculated according Reich et al, 1992 by :
!!               \latexonly
!!               \input{stomate_data_SLA.tex}
!!               \endlatexonly
!!               The sapling (young) biomass for trees and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_tree.tex}
!!               \endlatexonly
!!               The sapling biomass for grasses and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_grass.tex}
!!               \endlatexonly
!!               The critical stem diameter is given by the following formula :
!!               \latexonly
!!               \input{stomate_data_stem_diameter.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months)) 
!!                   \latexonly
!!                   \input{stomate_data_SLA_Reich_97.tex}
!!                   \endlatexonly
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!! - Reich PB, Walters MB, Ellsworth DS, (1992), Leaf life-span in relation to leaf, plant and 
!! stand characteristics among diverse ecosystems. Ecological Monographs, Vol 62, pp 365-392.
!! - Reich PB, Walters MB, Ellsworth DS (1997) From tropics to tundra: global convergence in plant 
!!  functioning. Proc Natl Acad Sci USA, 94:13730 13734
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE data


    INTEGER(i_std)                               :: i,j     !! Index (unitless)

!_ ================================================================================================================================

    ! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    printlev_loc=printlev
    !- pheno_gdd_crit
    pheno_gdd_crit(:,1) = pheno_gdd_crit_c(:)
    pheno_gdd_crit(:,2) = pheno_gdd_crit_b(:)         
    pheno_gdd_crit(:,3) = pheno_gdd_crit_a(:) 
    !
    !- senescence_temp
    senescence_temp(:,1) = senescence_temp_c(:)
    senescence_temp(:,2) = senescence_temp_b(:)
    senescence_temp(:,3) = senescence_temp_a(:)
    ! 
    !-LC 
    LC(:,ileaf) = LC_leaf(:) 
    LC(:,isapabove) = LC_sapabove(:) 
    LC(:,isapbelow) = LC_sapbelow(:) 
    LC(:,iheartabove) = LC_heartabove(:) 
    LC(:,iheartbelow) = LC_heartbelow(:) 
    LC(:,iroot) = LC_root(:) 
    LC(:,ifruit) = LC_fruit(:) 
    LC(:,icarbres) = LC_carbres(:) 
    LC(:,ilabile) = LC_labile(:) 

    IF ( printlev_loc >= 2 ) WRITE(numout,*) 'data: PFT characteristics'

    DO j = 2,nvm ! Loop over # PFTS 

       IF ( printlev_loc >= 2 ) WRITE(numout,'(a,i3,a,a)') '    > PFT#',j,': ', PFT_name(j)

       !
       ! 1 tree? (true/false)
       !
       IF ( printlev_loc >= 2 ) WRITE(numout,*) '       tree: (::is_tree) ', is_tree(j)


       !
       ! 4 specific leaf area per mass carbon = 2 * sla / dry mass (m^2.gC^{-1})
       !

       ! S. Zaehle: Reich et al, 1992 find no statistically significant differences between broadleaved and coniferous
       ! forests, specifically, the assumption that grasses grow needles is not justified. Replacing the function
       ! with the one based on Reich et al. 1997. Given that sla=100cm2/gDW at 9 months, sla is:
       ! sla=exp(5.615-0.46*ln(leaflon in months))

       ! Oct 2010 : sla values are prescribed by values given by N.Viovy 

       ! includes conversion from 
       !!       sla(j) = 2. * 1e-4 * EXP(5.615 - 0.46 * log(12./(longevity_eff_leaf(j)/365)))
       !!\latexonly
       !!\input{stomate_data_SLA.tex}
       !!\endlatexonly
!       IF ( leaf_tab(j) .EQ. 2 ) THEN
!
!          ! needle leaved tree
!          sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./(longevity_eff_leaf(j)/365)) ) ) *1e-4
!
!       ELSE
!
!          ! broad leaved tree or grass (Reich et al 1992)
!          sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./(longevity_eff_leaf(j)/365)) ) ) *1e-4
!
!       ENDIF

!!!$      IF ( leaf_tab(j) .EQ. 1 ) THEN
!!!$
!!!$        ! broad leaved tree
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./(longevity_eff_leaf(j)/365)) ) ) *1e-4
!!!$
!!!$      ELSE
!!!$
!!!$        ! needle leaved or grass (Reich et al 1992)
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./(longevity_eff_leaf(j)/365)) ) ) *1e-4
!!!$
!!!$      ENDIF
!!!$
!!!$      IF ( ( leaf_tab(j) .EQ. 2 ) .AND. ( pheno_type_tab(j) .EQ. 2 ) ) THEN
!!!$
!!!$        ! summergreen needle leaf
!!!$
!!!$        sla(j) = 1.25 * sla(j)
!!!$
!!!$      ENDIF

!!$       IF ( printlev_loc >= 2 ) WRITE(numout,*) '       specific leaf area (m**2/gC) (::sla):', sla(j), 12./(longevity_eff_leaf(j)/365)


       !+++CHECK+++
       ! May be needed for the DGVM but no longer needed for stomate_prescribe. There
       ! is a separate subroutine make_sapling that can be called whenever the initial
       ! biomass of a sapling needs to be calculated. Estimating sapling biomass 
       ! on-the-fly rather than once at the start of a simulation has the advantage that 
       ! allocation factors can be adapted to the site conditions.
       !+++++++++++

       !
       ! 6 migration speed (m/year)
       !

       IF ( is_tree(j) ) THEN

          migrate(j) = migrate_tree

       ELSE

          ! can be any value as grasses are, per *definition*, everywhere (big leaf).
          migrate(j) = migrate_grass

       ENDIF !( is_tree(j) )

       IF ( printlev_loc >= 2 ) WRITE(numout,*) '       migration speed (m/year): (::migrate(j))', migrate(j)

       !
       ! 7 critical stem diameter: beyond this diameter, the crown area no longer
       !     increases (m)
       !

       IF ( is_tree(j) ) THEN

          !!\latexonly
          !!\input{stomate_data_stem_diameter.tex}
          !!\endlatexonly          
          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year
         
       ELSE

          cn_sapl(j)=1

       ENDIF !( is_tree(j) )
   
       
       !
       ! 8 Coldest tolerable temperature (K)
       !

       IF ( ABS( tmin_crit(j) - undef ) .GT. min_stomate ) THEN
          tmin_crit(j) = tmin_crit(j) + ZeroCelsius
       ELSE
          tmin_crit(j) = undef
       ENDIF 

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       coldest tolerable temperature (K): (::tmin_crit(j))', tmin_crit(j)

       !
       ! 9 Maximum temperature of the coldest month: need to be below this temperature
       !      for a certain time to regrow leaves next spring *(vernalization)* (K)
       !

       IF ( ABS ( tcm_crit(j) - undef ) .GT. min_stomate ) THEN
          tcm_crit(j) = tcm_crit(j) + ZeroCelsius
       ELSE
          tcm_crit(j) = undef
       ENDIF

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       vernalization temperature (K): (::tcm_crit(j))', tcm_crit(j)

       !
       ! 10 critical values for phenology
       !

       ! 10.1 model used

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       phenology model used: (::pheno_model(j)) ',pheno_model(j)

       ! 10.2 growing degree days. What kind of gdd is meant (i.e. threshold 0 or -5 deg C
       !        or whatever), depends on how this is used in stomate_phenology.


       IF ( ( printlev_loc >= 2 ) .AND. ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical GDD is a function of long term T (C): (::gdd)'
          WRITE(numout,*) '          ',pheno_gdd_crit(j,1), &
               ' + T *',pheno_gdd_crit(j,2), &
               ' + T^2 *',pheno_gdd_crit(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( pheno_model(j) .EQ. 'moigdd' ) .OR. &
            ( pheno_model(j) .EQ. 'humgdd' )       ) .AND. &
            ( ANY(pheno_gdd_crit(j,:) .EQ. undef) )                      ) THEN
          CALL ipslerr_p(3,'stomate_data','problem with phenology parameters, critical GDD. (::pheno_model)','','')
       ENDIF

       ! 10.3 number of growing days

       IF ( ( printlev_loc >= 2 ) .AND. ( ngd_crit(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical NGD: (::ngd_crit(j))', ngd_crit(j)

       ! 10.4 critical temperature for ncd vs. gdd function in phenology (C)

       IF ( ( printlev_loc >= 2 ) .AND. ( ncdgdd_temp(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical temperature for NCD vs. GDD (C): (::ncdgdd_temp(j))', &
            ncdgdd_temp(j)

       ! 10.5 humidity fractions (0-1, unitless)

       IF ( ( printlev_loc >= 2 ) .AND. ( hum_frac(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical humidity fraction: (::hum_frac(j))', &
            &  hum_frac(j)


       ! 10.6 minimum time elapsed since moisture minimum (days)

       IF ( ( printlev_loc >= 2 ) .AND. ( hum_min_time(j) .NE. undef ) ) &
            WRITE(numout,*) '         time to wait after moisture min (d): (::hum_min_time(j))', &
        &    hum_min_time(j)

       !
       ! 11 critical values for senescence
       !

       ! 11.1 type of senescence

       IF ( printlev_loc >= 2 ) WRITE(numout,*) '       type of senescence: (::senescence_type(j))',&
            senescence_type(j)

       ! 11.2 critical temperature for senescence (C)

       IF ( ( printlev_loc >= 2 ) .AND. ( ALL(senescence_temp(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical temperature for senescence (C) is'
          WRITE(numout,*) '          a function of long term T (C): (::senescence_temp)'
          WRITE(numout,*) '          ',senescence_temp(j,1), &
               ' + T *',senescence_temp(j,2), &
               ' + T^2 *',senescence_temp(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'cold' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )      ) .AND. &
            ( ANY(senescence_temp(j,:) .EQ. undef ) )           ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, temperature. (::senescence_type)','','')
       ENDIF

       ! 11.3 critical relative moisture availability for senescence

       IF ( ( printlev_loc >= 2 ) .AND. ( senescence_hum(j) .NE. undef ) ) THEN
          WRITE(numout,*)  ' max. critical relative moisture availability for' 
          WRITE(numout,*)  ' senescence: (::senescence_hum(j))',  &
               & senescence_hum(j)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( senescence_hum(j) .EQ. undef )                   ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, humidity.(::senescence_type)','','')
       ENDIF

       ! 14.3 relative moisture availability above which there is no moisture-related
       !      senescence (0-1, unitless)

       IF ( ( printlev_loc >= 2 ) .AND. ( nosenescence_hum(j) .NE. undef ) ) THEN
          WRITE(numout,*) '         relative moisture availability above which there is' 
          WRITE(numout,*) '             no moisture-related senescence: (::nosenescence_hum(j))', &
               &  nosenescence_hum(j)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( nosenescence_hum(j) .EQ. undef )                   ) THEN
          CALL ipslerr_p(3,'stomate_data','Problem with senescence parameters, humidity. (::senescence_type)','','')
       ENDIF

       !
       ! 12 sapwood -> heartwood conversion time (days)
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       sapwood -> heartwood conversion time (d): (::longevity_sap(j))', longevity_sap(j)

       !
       ! 13 fruit lifetime (days)
       !

       IF ( printlev_loc >= 2 ) WRITE(numout,*) '       fruit lifetime (d): (::longevity_fruit(j))', longevity_fruit(j)

       !
       ! 14 length of leaf death (days)
       !      For evergreen trees, this variable determines the lifetime of the leaves.
       !      Note that it is different from the value given in (longevity_eff_leaf/365).
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       length of leaf death (d): (::leaffall(j))', leaffall(j)


       !
       ! 17 minimum lai, initial (m^2.m^{-2})
       !

       IF ( is_tree(j) ) THEN
          lai_initmin(j) = lai_initmin_tree
       ELSE
          lai_initmin(j) = lai_initmin_grass
       ENDIF !( is_tree(j) )

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       initial LAI: (::lai_initmin(j))', lai_initmin(j)

       !
       ! 19 maximum LAI (m^2.m^{-2})
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       critical LAI above which no leaf allocation: (::lai_max(j))', lai_max(j)

       !
       ! 20 fraction of primary leaf and root allocation put into reserve (0-1, unitless)
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       reserve allocation factor: (::ecureuil(j))', ecureuil(j)


       ! 23 natural ?
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       Natural: (::natural(j))', natural(j)

       !
       ! 24 Vcmax et Vjmax (umol.m^{-2}.s^{-1}) 
       !

       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       Maximum rate of carboxylation: (::Vcmax_25(j))', vcmax25(j)
       !
       ! 25 constants for photosynthesis temperatures
       !

       IF ( printlev_loc >= 2 ) THEN


          !
          ! 26 Properties
          !
          WRITE(numout,*) '       C4 photosynthesis: (::is_c4(j))', is_c4(j)

       ENDIF

       !
       ! 27 extinction coefficient of the Monsi and Saeki (1953) relationship 
       !
       IF ( printlev_loc >= 2 ) THEN
          WRITE(numout,*) '       extinction coefficient: (::ext_coeff(j))', ext_coeff(j)
       ENDIF

       !
       ! 30 fraction of allocatable biomass which is lost as growth respiration (0-1, unitless)
       !
       IF ( printlev_loc >= 2 ) &
            WRITE(numout,*) '       growth respiration fraction: (::frac_growthresp(j))', frac_growthresp(j)

    ENDDO ! Loop over # PFTS 

    !
    ! 29 time scales for phenology and other processes (in days)
    !

    tau_longterm_max = coeff_tau_longterm * one_year

    IF ( printlev_loc >= 2 ) THEN

       WRITE(numout,*) '   > time scale for ''monthly'' moisture availability (d): (::tau_hum_month)', &
            tau_hum_month
       WRITE(numout,*) '   > time scale for ''weekly'' moisture availability (d): (::tau_hum_week)', &
           tau_hum_week
       WRITE(numout,*) '   > time scale for ''monthly'' 2 meter temperature (d): (::tau_t2m_month)', &
            tau_t2m_month
       WRITE(numout,*) '   > time scale for ''weekly'' 2 meter temperature (d): (::tau_t2m_week)', &
            tau_t2m_week
       WRITE(numout,*) '   > time scale for ''weekly'' GPP (d): (::tau_gpp_week)', &
            tau_gpp_week
       WRITE(numout,*) '   > time scale for ''monthly'' soil temperature (d): (::tau_tsoil_month)', &
            tau_tsoil_month
       WRITE(numout,*) '   > time scale for vigour calculations (y): (::tau_longterm_max / one_year)', &
            tau_longterm_max / one_year

    ENDIF

    !
    ! 30 Maintenance respiration
    !

    maint_resp_slope(:,1) = maint_resp_slope_c(:)               
    maint_resp_slope(:,2) = maint_resp_slope_b(:) 
    maint_resp_slope(:,3) = maint_resp_slope_a(:)
 
    IF (printlev_loc >= 4) WRITE(numout,*) 'Leaving stomate_data'

  END SUBROUTINE data


!! ================================================================================================================================
!! SUBROUTINE 	: calculate_aed
!!
!>\BRIEF         Compute the Average Edge Distance (AED) used by SPITFIRE and the canopy
!!               light edge-effect (LIGHT) submodules.
!!
!! DESCRIPTION : Centralises the fragmentation geometry. Two AED variables coexist :
!!               AED_fire (denominator = area_burnable_veg, used by stomate_spitfire) and
!!               AED_light (denominator = area_forest = is_tree only, used by stomate_laieff).
!!               edge_length (m) is the shared input, read annually from Edge_length.nc.
!!
!! REFERENCE(S) : Bowring et al. (2024) ; Marie (2026) edge_effect.pdf.
!! \n
!_ ================================================================================================================================

  SUBROUTINE calculate_aed(npts, veget_max, area, contfrac, edge_length)

    !! 0. Arguments
    INTEGER(i_std), INTENT(in)                    :: npts          !! Number of grid points
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  :: veget_max     !! Max vegetation fraction per PFT (unitless)
    REAL(r_std), DIMENSION(npts),     INTENT(in)  :: area          !! Grid cell area (m^2)
    REAL(r_std), DIMENSION(npts),     INTENT(in)  :: contfrac      !! Continental fraction (unitless)
    REAL(r_std), DIMENSION(npts),     INTENT(in)  :: edge_length   !! Edge length from Edge_length.nc (m)

    !! 0.2 Local variables
    INTEGER(i_std)                                :: ivm, ier
    REAL(r_std), DIMENSION(npts)                  :: area_burnable_veg !! ha
    REAL(r_std), DIMENSION(npts)                  :: area_land         !! ha
    REAL(r_std), DIMENSION(npts)                  :: burnable_vegfrac, forest_vegfrac

    !! 1. Early return if no consumer of AED_fire/AED_light is active.
    !! Guillaume M. -- ok_spitfire MUST be part of the test: SPITFIRE uses AED_fire
    !! unconditionally (fedge_wind, calculate_human_ignition, combustion_fraction),
    !! so a plain fire run with every OK_AED_* off would leave it unallocated.
    IF (.NOT. (ok_spitfire    .OR. &
               ok_aed_humign  .OR. ok_aed_fuel .OR. ok_aed_size .OR. &
               ok_aed_wind    .OR. ok_aed_light)) RETURN

    !! 2. Lazy allocation of module-level outputs
    IF (.NOT. ALLOCATED(AED_fire)) THEN
       ALLOCATE(AED_fire(npts), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'calculate_aed','Pb alloc AED_fire','','')
       AED_fire(:) = 50000.0
    ENDIF
    IF (.NOT. ALLOCATED(AED_light)) THEN
       ALLOCATE(AED_light(npts), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'calculate_aed','Pb alloc AED_light','','')
       AED_light(:) = 50000.0
    ENDIF
    IF (.NOT. ALLOCATED(area_forest)) THEN
       ALLOCATE(area_forest(npts), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'calculate_aed','Pb alloc area_forest','','')
       area_forest(:) = zero
    ENDIF

    !! 3. Land area in hectares (matches stomate_spitfire.f90:471)
    area_land(:) = area(:) * contfrac(:) / 10000.0

    !! 4. Burnable fraction (trees + grass that can burn) — SPITFIRE denominator
    burnable_vegfrac(:) = zero
    DO ivm = 1, nvm
       IF (burnable(ivm)) burnable_vegfrac(:) = burnable_vegfrac(:) + veget_max(:,ivm)
    ENDDO
    area_burnable_veg(:) = area_land(:) * burnable_vegfrac(:)

    !! 5. Forest fraction (is_tree only) — LIGHT denominator
    forest_vegfrac(:) = zero
    DO ivm = 1, nvm
       IF (is_tree(ivm)) forest_vegfrac(:) = forest_vegfrac(:) + veget_max(:,ivm)
    ENDDO
    area_forest(:) = area_land(:) * forest_vegfrac(:)

    !! 6. AED_fire — preserves the historical SPITFIRE behaviour
    WHERE (edge_length(:) > min_stomate)
       AED_fire(:) = area_burnable_veg(:) * 10000.0 / edge_length(:) / 2.0
    ELSEWHERE
       AED_fire(:) = 50000.0
    ENDWHERE
    AED_fire(:) = MIN(50000.0_r_std, AED_fire(:))
    WHERE (area_burnable_veg(:) < min_stomate) AED_fire(:) = 50000.0

    IF ((ok_aed_humign .OR. ok_aed_fuel .OR. ok_aed_size .OR. ok_aed_wind) &
         .AND. ANY(AED_fire(:) < edge_distance_wetness/2)) THEN
       WRITE(numout,*) 'edge_length = ', edge_length(:)
       WRITE(numout,*) 'area_burnable_veg = ', area_burnable_veg(:)
       WRITE(numout,*) 'edge_distance_wetness/2 = ', edge_distance_wetness/2
       WRITE(numout,*) 'AED_fire = ', AED_fire(:)
       CALL ipslerr_p(3, 'calculate_aed', 'AED_fire < edge_distance_wetness/2', '', '')
    ENDIF

    !! 7. AED_light — new denominator for canopy light edge effect
    WHERE (edge_length(:) > min_stomate)
       AED_light(:) = area_forest(:) * 10000.0 / edge_length(:) / 2.0
    ELSEWHERE
       AED_light(:) = 50000.0
    ENDWHERE
    AED_light(:) = MIN(50000.0_r_std, AED_light(:))
    WHERE (area_forest(:) < min_stomate) AED_light(:) = 50000.0

    ! Guillaume M. -- Lower clamp at edge_distance_light/2. Below it the edge penetration
    ! covers the whole patch and fedge_pgap saturates at 1. Clamping is the physically
    ! consistent answer, so it is done silently: a fatal error here fires on every global
    ! run holding cells of small forest fraction. A level-2 warning flags the setup.
    IF (ok_aed_light) THEN
       IF (ANY(AED_light(:) < edge_distance_light/2)) THEN
          CALL ipslerr_p(2, 'calculate_aed',                                       &
               'AED_light < edge_distance_light/2 on some grid cells —',          &
               'clamping at the lower threshold (fedge_pgap saturates).',         &
               '')
       ENDIF
       AED_light(:) = MAX(edge_distance_light/2, AED_light(:))
    ENDIF

  END SUBROUTINE calculate_aed


!! ================================================================================================================================
!! SUBROUTINE    : update_edge_length
!!
!>\BRIEF          Advance edge_length_dyn by one timestep from perturbation
!!                accumulators dA_fire/dA_storm/dA_pest, with exponential
!!                relaxation toward the continuous-forest baseline EL_base =
!!                aed_edge_ref_frac*edge_length_ref (default frac=0).
!!
!! DESCRIPTION   : Feedback baseline EL_base = aed_edge_ref_frac*edge_length_ref
!!                 avoids double-counting the EFDA observation (which already holds
!!                 the historical-perturbation edge); frac=0 -> pure perturbation-
!!                 driven edge. Two paths. (a) Legacy AED_FEEDBACK (ok_edge_from_age_class=
!!                 FALSE): forward Euler d(EL)/dt = Sum_i alpha_i*dA_i/dt - (EL-EL_base)/tau_rec.
!!                 (b) AED_REGROWTH (ok_edge_from_age_class=TRUE): the dA_* are ingested
!!                 into a class-0 conveyor of annual cohorts that graduate after
!!                 n_recruit years; edge_length_dyn = EL_base + Sum_i alpha_i *
!!                 class0_area_i (mechanistic recovery, no tau_rec). veget_max is
!!                 NOT modified (Level-1). Either way edge_length_dyn is bounded by
!!                 the multi-patch cap and the dA_* are reset on exit.
!!                 No-op when ok_aed_feedback is FALSE.
!_ ================================================================================================================================

  SUBROUTINE set_management_intensity(npts, veget_max)
    !! Guillaume M. -- Single entry point of forest management
    !! (MODULE_DESIGN_MANAGEMENT_INTENSITY.md). From the intensity-class fractions mi_frac
    !! alone it derives target_rotation_age (the ONLY channel by which intensity acts on
    !! the dynamics), alpha_harvest (harvest edge coefficient) and clearcut_area (realised
    !! patch size, edge cap and diagnostic only). Called once per step, before the edge update.
    INTEGER(i_std),                      INTENT(in) :: npts       !! Nombre de points de grille
    REAL(r_std),    DIMENSION(npts,nvm), INTENT(in) :: veget_max  !! Fraction max de vegetation par PFT

    ! --- Variables locales
    INTEGER(i_std)                    :: ji, jv, ic, icls
    REAL(r_std)                       :: ftot, s_wmean, a_wmean, hmod, d_wmean
    REAL(r_std)                       :: acls_tot, jmat
    REAL(r_std), DIMENSION(nmiclass)  :: alpha_c    !! alpha par classe = k/sqrt(S_c) (m/m2)
    REAL(r_std), DIMENSION(nmiclass)  :: fc         !! fractions de classe de la maille, repli applique (-)
    REAL(r_std), DIMENSION(nagec)     :: acls       !! surface par classe d'age dans la maille
    REAL(r_std)                       :: r_wmean    !! moyenne ponderee de 1/mi_rotation_factor (-)
    REAL(r_std), DIMENSION(nvm)       :: rref       !! reference rotation of this cell: map where covered, namelist elsewhere (yr)

    IF (.NOT. ok_management_intensity) RETURN
    IF (.NOT. ALLOCATED(mi_frac)) RETURN

    IF (.NOT. ALLOCATED(clearcut_area))          ALLOCATE(clearcut_area(npts))
    IF (.NOT. ALLOCATED(alpha_harvest))      ALLOCATE(alpha_harvest(npts))
    IF (.NOT. ALLOCATED(target_rotation_age)) ALLOCATE(target_rotation_age(npts,nvm))

    !! Guillaume M. -- Per-class alpha, spatially invariant, so computed once per call.
    !! alpha is then averaged LINEARLY over classes rather than derived from the averaged
    !! patch size: alpha ~ 1/sqrt(S) is convex, so averaging S first underestimates it.
    DO ic = 1, nmiclass
       alpha_c(ic) = mi_edge_shape_coef / SQRT(MAX(mi_clearcut_size(ic), min_stomate))
    END DO

    DO ji = 1, npts

       !! Guillaume M. -- Regional reference where ROTATION_REF_FILE has one, namelist scalar
       !! elsewhere (design/MODULE_DESIGN_ROTATION_REGION.md). The map replaces the REFERENCE
       !! only: the intensity factor still applies below. Not allocated = NONE = bit-neutral.
       rref(:) = rotation_ref(:)
       IF (ALLOCATED(rotation_ref_map)) THEN
          WHERE (rotation_ref_map(ji,:) > zero) rref(:) = rotation_ref_map(ji,:)
       ENDIF

       !! --- 1. Guillaume M. -- Normalise the class fractions: their sum can differ from 1
       !!        after interpolation, or on a cell holding no managed forest.
       ftot = zero
       DO ic = 1, nmiclass
          fc(ic) = MAX(mi_frac(ji,ic), zero)
          ftot   = ftot + fc(ic)
       END DO
       !! Guillaume M. -- Silent map: the gap is the edge of the EFDA footprint, where the
       !! cells still carry managed forest, so fall back to MI_DEFAULT_CLASS instead of
       !! claiming no harvest. Routed through the normal path, so the rotation, the edge
       !! coefficient, the patch size and the diameter stay consistent with each other.
       IF (ftot <= min_stomate) THEN
          fc(:)                 = zero
          fc(mi_default_class) = un
          ftot                  = un
       END IF

       !! --- 2. Guillaume M. -- Means weighted by the managed area.
       a_wmean = zero
       s_wmean = zero
       DO ic = 1, nmiclass
          a_wmean = a_wmean + fc(ic) * alpha_c(ic)          ! moyenne LINEAIRE de alpha
          s_wmean = s_wmean + fc(ic) * mi_clearcut_size(ic)
       END DO
       !!
       !! Guillaume M. -- Harvest cadence is PER PFT: the rate is an absolute flux, so a
       !! class-only rate would turn an oak stand (100-150 yr) and a eucalypt (10-35 yr) at
       !! the same speed. /!\ The additive quantity is the RATE, not the rotation.
       !!   rate(ji,jv) = (1/rref(jv)) * Sigma_c f_c/mi_rotation_factor(c) / ftot
       !!   with rref = ROTATION_REF(jv), or the regional map where ROTATION_REF_FILE covers the cell.
       !! The class factor is PFT independent: averaged ONCE, then divided per PFT.
       r_wmean = zero
       DO ic = 1, nmiclass
          r_wmean = r_wmean + fc(ic) / mi_rotation_factor(ic)
       END DO
       r_wmean = r_wmean / ftot
       !! Guillaume M. -- Zero means "no rotation defined": only tree PFTs read this map,
       !! and every tree has rotation_ref > 0, so they all take the WHERE branch below.
       target_rotation_age(ji,:) = zero
       WHERE (rref(:) > zero .AND. r_wmean > min_stomate) &
            target_rotation_age(ji,:) = rref(:) / r_wmean
       !! Guillaume M. -- ROTATION_GROWTH: blend the published itinerary with the rotation
       !! the stand can actually deliver. w = 0 keeps the recommendation alone (bit-neutral
       !! default), w = 1 follows the site. A CONVEX combination, so the result can never
       !! leave the interval spanned by the two terms. rotation_growth <= 0 means "not
       !! estimable here" -- the recommendation is then used untouched, never extrapolated.
       IF (rotation_growth_weight > zero .AND. ALLOCATED(rotation_growth)) THEN
          WHERE (target_rotation_age(ji,:) > zero .AND. rotation_growth(ji,:) > zero) &
               target_rotation_age(ji,:) = &
                    (un - rotation_growth_weight) * target_rotation_age(ji,:) &
                    + rotation_growth_weight * rotation_growth(ji,:)
       ENDIF
       alpha_harvest(ji)       = a_wmean / ftot
       s_wmean                     = s_wmean / ftot

       !! Guillaume M. -- Clearcut diameter modulated by intensity: weighted mean of the
       !! class factors. It is a FACTOR and not an absolute diameter, so the difference
       !! between species is kept while the target moves with management.
       IF (ALLOCATED(dia_factor)) THEN
          d_wmean = zero
          DO ic = 1, nmiclass
             d_wmean = d_wmean + fc(ic) * mi_dia_factor(ic)
          END DO
          dia_factor(ji) = d_wmean / ftot
       END IF

       !! --- 3. Guillaume M. -- Maturity J of the age structure (0 = all young, 1 = all
       !!        mature). Where mature forest forms a large block the realised cuts are
       !!        large, where it is young or scattered they are small. J is directional,
       !!        unlike a concentration index, so clearcut_area becomes a PHASE diagnostic.
       acls(:)  = zero
       acls_tot = zero
       DO jv = 1, nvm
          IF (.NOT. is_tree(jv)) CYCLE
          icls = jv - start_index(agec_group(jv)) + 1        ! rang de classe d'age 1..nagec
          IF (icls >= 1 .AND. icls <= nagec) THEN
             acls(icls) = acls(icls) + veget_max(ji,jv)
             acls_tot   = acls_tot   + veget_max(ji,jv)
          END IF
       END DO

       hmod = un
       IF (acls_tot > min_stomate .AND. nagec > 1) THEN
          jmat = zero
          DO icls = 1, nagec
             jmat = jmat + acls(icls) * REAL(icls - 1, r_std) / REAL(nagec - 1, r_std)
          END DO
          jmat = jmat / acls_tot
          !! h(J) = 1 + strength * (J - J_ref)/(1 - J_ref), J_ref = 0.5 (balanced structure)
          hmod = un + mi_age_mod_strength * (jmat - 0.5_r_std) / 0.5_r_std
          hmod = MAX(hmod, zero)
       END IF

       !! --- 4. Realised patch size, bounded.
       clearcut_area(ji) = MIN(MAX(s_wmean * hmod, mi_clearcut_size_min), mi_clearcut_size_max)

    END DO

  END SUBROUTINE set_management_intensity


  SUBROUTINE update_edge_length(npts, dt, area, contfrac, veget_max)

    !! 0. Arguments
    INTEGER(i_std), INTENT(in)                    :: npts      !! Number of grid points
    REAL(r_std),    INTENT(in)                    :: dt        !! Timestep duration (s)
    REAL(r_std), DIMENSION(npts),     INTENT(in)  :: area      !! Grid cell area (m²)
    REAL(r_std), DIMENSION(npts),     INTENT(in)  :: contfrac  !! Continental fraction (-)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  :: veget_max !! Max vegetation fraction per PFT (AED_REGROWTH: per-PFT recovery weighting)

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(npts)                  :: dEL_source, dEL_recovery, max_edge_phys
    REAL(r_std), DIMENSION(npts)                  :: alpha_harvest_use  !! coef de lisiere harvest par maille (intensite de gestion, m/m2)
    REAL(r_std), DIMENSION(npts)                  :: patch_area_eff     !! taille de parcelle pour le cap de lisiere (m2)
    REAL(r_std), DIMENSION(npts)                  :: alpha_fire_use     !! coef de lisiere feu par maille (m/m2)
    REAL(r_std), DIMENSION(npts)                  :: alpha_storm_use    !! coef de lisiere tempete par maille (m/m2)
    REAL(r_std), DIMENSION(npts)                  :: alpha_pest_use     !! coef de lisiere scolyte par maille (m/m2)
    REAL(r_std), DIMENSION(npts)                  :: edge_src, class0_tot, frac_class0_diag
    REAL(r_std), DIMENSION(npts)                  :: edge_length_base   !! Continuous-forest baseline of the feedback (m)
    INTEGER(i_std)                                :: ib
    INTEGER(i_std)                                :: ji, jv, kk         !! AED_REGROWTH Trecruit(PFT) loop indices
    REAL(r_std)                                   :: w_for, w_rec       !! per-point forest veget sum and recovery-weighted sum
    REAL(r_std), DIMENSION(npts)                  :: A_open, alpha_eff, flux_tot   !! AED_REGROWTH v2 (A_open from young age class)
    REAL(r_std), DIMENSION(npts)                  :: f_class1                      !! Age-class-1 cover fraction, openness-weighted if asked (-)
    REAL(r_std), DIMENSION(npts,nvm)              :: rref_out                      !! ROTATION_REF_MAP diagnostic: reference rotation used (yr)
    REAL(r_std), DIMENSION(npts,nvm)              :: mtt_out                       !! MAT_TARGET_TERM diagnostic: effective terminal setpoint (-)
    REAL(r_std), DIMENSION(npts,nvm)              :: mae_out                       !! MAT_AGE_ENTRY / ROTATION_GROWTH diagnostics (yr)
    REAL(r_std), DIMENSION(npts)                  :: edge_meso, grain_eff          !! FRAG_3TERMS: meso edge (m) and per-cell grain (m)
    REAL(r_std), DIMENSION(npts)                  :: edge_macro                    !! FRAG_3TERMS: macro edge from roads (m)
    REAL(r_std), DIMENSION(npts)                  :: edge_hydro                    !! FRAG_3TERMS: permanent edge from lakes, coast and rivers (m)
    REAL(r_std), DIMENSION(npts)                  :: edge_micro                    !! FRAG_3TERMS: micro edge, what the dynamics add above the prescribed base (m)
    REAL(r_std)                                   :: w_grain                       !! FRAG_3TERMS: veget_max-weighted grain accumulator
    REAL(r_std), DIMENSION(npts)                  :: alpha_bg                      !! v2 piste B: climatological background α per cell
    REAL(r_std)                                   :: decay_mix                     !! v2 agent-flux decay
    INTEGER(i_std)                                :: icls                          !! v2 age-class index of a PFT
    !! Guillaume M. -- Climatological agent mix used for the background edge of young or
    !! open forest carrying no current disturbance flux (weights from EFDA attribution).
    REAL(r_std), PARAMETER :: bgw_fire=0.045_r_std, bgw_storm=0.07_r_std, bgw_pest=0.02_r_std, bgw_harvest=0.86_r_std

    IF (.NOT. ok_aed_feedback) RETURN
    IF (.NOT. (ALLOCATED(edge_length_dyn) .AND. ALLOCATED(edge_length_ref))) RETURN

    !! 0.3 Guillaume M. -- Feedback baseline. edge_length_ref is the OBSERVED EFDA edge and
    !!     already integrates the historical perturbations, so using it as additive floor or
    !!     relaxation target double-counts that history. The anchor is instead a continuous
    !!     forest baseline = aed_edge_ref_frac·edge_length_ref: frac=0 builds EL_dyn purely
    !!     from the simulated source, frac=1 restores the EFDA-anchored behaviour.
    edge_length_base(:) = aed_edge_ref_frac * edge_length_ref(:)

    !! 0bis. Guillaume M. -- MESO-SCALE term (MODULE_DESIGN_FRAGMENTATION_3TERMES.md).
    !!   EL_meso = A_cell * 4*f*(1-f) / r, f the forest cover fraction, r the sub-grid
    !!   grain: zero edge at f=0 and f=1, maximum at f=0.5. This is the edge a mature
    !!   contiguous forest carries against the non-forest matrix, absent from the micro
    !!   term which only sees disturbance. r is per-MTC and needs no map.
    edge_meso(:) = zero
    DO ji = 1, npts
       w_for = zero; w_grain = zero
       DO jv = 2, nvm
          IF (.NOT. is_tree(jv)) CYCLE
          w_for   = w_for   + veget_max(ji,jv)
          w_grain = w_grain + veget_max(ji,jv) * edge_grain_mtc_of(jv)
       ENDDO
       IF (w_for > min_stomate) THEN
          grain_eff(ji) = MAX(w_grain / w_for, min_stomate)
          edge_meso(ji) = area(ji) * contfrac(ji) * quatre * w_for * (un - w_for) &
                          / grain_eff(ji)
       ELSE
          grain_eff(ji) = edge_grain_default
       ENDIF
    ENDDO
    edge_length_base(:) = edge_length_base(:) + edge_meso(:)

    !! 0ter. Guillaume M. -- MACRO-SCALE term, roads. EL_macro = 2 * RL_forest * c_GRIP.
    !!   The 2 is GEOMETRIC: a road of length L cuts 2L of forest edge, one on each side.
    !!   /!\ It is NOT the other factor 2 of Bowring et al. (2024), which corrects GRIP
    !!   under-reporting and lives in road_grip_correction; conflating them costs a factor 2.
    !!   Additive with the meso term, which does not see roads: at 30 m they hide under canopy.
    edge_macro(:) = zero
    IF (ALLOCATED(road_length_map)) THEN
       IF (ok_road_forest_mask) THEN
          ! The map already carries the 5 arcmin forest restriction (EFDA 30 m masks).
          edge_macro(:) = deux * road_length_map(:) * road_grip_correction
       ELSE
          ! Raw road length: restrict with the model's own forest fraction, so the term
          ! follows land-use change -- at the price of overestimating, roads being denser
          ! outside forest than inside.
          DO ji = 1, npts
             w_for = zero
             DO jv = 2, nvm
                IF (is_tree(jv)) w_for = w_for + veget_max(ji,jv)
             ENDDO
             edge_macro(ji) = deux * road_length_map(ji) * road_grip_correction &
                              * MIN(un, MAX(zero, w_for))
          ENDDO
       ENDIF
    ENDIF
    edge_length_base(:) = edge_length_base(:) + edge_macro(:)

    !! 0quater. Guillaume M. -- PERMANENT term, water. Lakes, coastline and rivers cut the
    !!   forest exactly as roads do, and the meso term cannot see them: it is built from
    !!   veget_max on area_land, whose non-forest side is TERRESTRIAL only.
    !!   /!\ Added AS IS. The map already carries its geometric factors -- the two banks of a
    !!   river are counted in it -- and is already restricted to forest. Re-applying a factor 2
    !!   or a forest weighting here would double-count, as for the roads above.
    edge_hydro(:) = zero
    IF (ALLOCATED(edge_hydro_map)) edge_hydro(:) = edge_hydro_map(:)
    edge_length_base(:) = edge_length_base(:) + edge_hydro(:)


    !! 1. Guillaume M. -- Per-agent edge-production coefficients alpha (m/m2). Harvest alpha
    !!    follows management intensity, hence it is spatial and comes from
    !!    set_management_intensity. Fire, storm and pest are near-uniform in space, so they
    !!    are single calibrated values. /!\ mi_edge_shape_coef is tuned against the forest
    !!    area-weighted harvest alpha: changing one without the other breaks the constraint.
    IF (ok_management_intensity .AND. ALLOCATED(alpha_harvest)) THEN
       alpha_harvest_use(:) = alpha_harvest(:)
       patch_area_eff(:)    = MAX(clearcut_area(:), min_stomate)
    ELSE
       alpha_harvest_use(:) = alpha_harvest_edge
       patch_area_eff(:)    = edge_min_patch_area
    END IF

    alpha_fire_use(:)  = alpha_fire_edge
    alpha_storm_use(:) = alpha_storm_edge
    alpha_pest_use(:)  = alpha_pest_edge

    IF (ok_edge_from_age_class .AND. ALLOCATED(class0_area)) THEN
       !! 1a-3a. Guillaume M. -- AED_REGROWTH path, class-0 conveyor. Ingest this step's
       !!        stand-replaced area into the youngest cohort, advance the conveyor once a
       !!        year, then DERIVE edge_length_dyn from the live stock (mechanistic
       !!        recovery, no tau_rec). veget_max is NOT modified here.

       IF (ALLOCATED(agent_flux)) THEN
          !! === AED_REGROWTH V2 (MODULE_DESIGN_AED_REGROWTH_V2.md) ===
          !! Guillaume M. -- Edge from the YOUNG age-class area A_open, which the NAGEC
          !! machinery keeps below the forest area, so FRAC_CLASS0 <= 1 and nothing is
          !! counted twice. Recovery time is the age-class residence, PFT dependent through
          !! growth. Agent mix from a leaky dA-flux integrator of memory aed_mix_tau.
          decay_mix = EXP( - dt / MAX(aed_mix_tau * one_year * one_day, dt) )
          agent_flux(:,IDX_FIRE)    = agent_flux(:,IDX_FIRE)    * decay_mix + dA_fire(:)
          agent_flux(:,IDX_STORM)   = agent_flux(:,IDX_STORM)   * decay_mix + dA_storm(:)
          agent_flux(:,IDX_PEST)    = agent_flux(:,IDX_PEST)    * decay_mix + dA_pest(:)
          agent_flux(:,IDX_HARVEST) = agent_flux(:,IDX_HARVEST) * decay_mix + dA_harvest(:)

          !! Guillaume M. -- A_open is the regeneration gap and nothing else: age class 1,
          !! the sole genuine post-disturbance opening. This is a membership test, not a
          !! weighting. The forest/non-forest boundary the mature classes once carried is
          !! produced by the meso term of edge_length_base from the landscape geometry --
          !! a second f*(1-f) term here would count the same interface twice.
          f_class1(:) = zero
          DO jv = 1, nvm
             IF (.NOT. is_tree(jv)) CYCLE
             icls = jv - start_index(agec_group(jv)) + 1   ! age class 1..nagec
             IF (icls /= 1) CYCLE
             IF (ok_aed_edge_rdi_weight .AND. ALLOCATED(rdi_stand)) THEN
                !! Guillaume M. -- Weight the class-1 area by the actual openness of the
                !! stand: without it a closing stand still counts as fully open. The weight
                !! (1-RDI) is high just after a cut and decays as the stand stocks up, so
                !! canopy reclosure becomes EMERGENT and subsumes tau_rec.
                !! /!\ RDI exceeds 1 on dense stands: bound it or it removes edge.
                f_class1(:) = f_class1(:) + veget_max(:,jv) &
                     * MAX(zero, un - MIN(un, rdi_stand(:,jv)))
             ELSE
                f_class1(:) = f_class1(:) + veget_max(:,jv)
             ENDIF
          ENDDO
          A_open(:) = aed_alpha_reg * f_class1(:) * area(:) * contfrac(:)   ! fraction -> m²

          ! agent-weighted α (per cell) from the flux mix: Σ_a α_a·flux_a / Σ flux
          flux_tot(:) = SUM(agent_flux, DIM=2)
          alpha_eff(:) =   alpha_fire_use(:)    * agent_flux(:,IDX_FIRE)    &
                         + alpha_storm_use(:)   * agent_flux(:,IDX_STORM)   &
                         + alpha_pest_use(:)    * agent_flux(:,IDX_PEST)    &
                         + alpha_harvest_use(:) * agent_flux(:,IDX_HARVEST)
          !! Guillaume M. -- Background alpha (MODULE_DESIGN_AED_REGROWTH_V2.md): a cell with
          !! young or open forest but no current disturbance flux still carries a standing
          !! fragmentation edge that the observations see. When flux_tot is null the agent
          !! mix is undefined, so A_open is converted with a fixed climatological mix rather
          !! than left at zero edge.
          alpha_bg(:) =   bgw_fire    * alpha_fire_use(:)    &
                        + bgw_storm   * alpha_storm_use(:)   &
                        + bgw_pest    * alpha_pest_use(:)    &
                        + bgw_harvest * alpha_harvest_use(:)
          WHERE (flux_tot(:) > min_stomate)
             edge_src(:) = A_open(:) * alpha_eff(:) / flux_tot(:)
          ELSEWHERE
             edge_src(:) = aed_bg_frac * A_open(:) * alpha_bg(:)   ! piste B: standing-fragmentation background edge
          ENDWHERE
          class0_tot(:) = A_open(:)   ! FRAC_CLASS0 diagnostic = A_open/forest (<= 1)

       ENDIF   ! agent_flux allocated
       edge_length_dyn(:) = edge_length_base(:) + edge_src(:)

    ELSE
       !! 1b-3b. Legacy AED_FEEDBACK path — additive source + tau_rec relaxation.
       dEL_source(:) =   alpha_fire_use(:)    * dA_fire(:)    &
                       + alpha_storm_use(:)   * dA_storm(:)   &
                       + alpha_pest_use(:)    * dA_pest(:)    &
                       + alpha_harvest_use(:) * dA_harvest(:)

       IF (tau_rec_edge > min_stomate) THEN
          dEL_recovery(:) = (edge_length_dyn(:) - edge_length_base(:)) * dt / tau_rec_edge
       ELSE
          dEL_recovery(:) = zero
       ENDIF

       edge_length_dyn(:) = edge_length_dyn(:) + dEL_source(:) - dEL_recovery(:)
    ENDIF

    !! 4. Guillaume M. -- Physical bound: maximum edge of a cell fully tiled by square
    !!    patches of area edge_min_patch_area, 4·(area·contfrac)/sqrt(A_min). A single-patch
    !!    cap 4·sqrt(area) under-bounds a fragmented landscape by 2-3 orders of magnitude,
    !!    a 1° cell holding of the order of 1e6 sub-hectare patches.
    max_edge_phys(:) = 4.0_r_std * MAX(zero, area(:) * contfrac(:)) &
                       / SQRT(MAX(patch_area_eff(:), min_stomate))
    edge_length_dyn(:) = MIN(max_edge_phys(:), MAX(zero, edge_length_dyn(:)))

    !! 4.5 Guillaume M. -- There is no automatic AED_SPINUP convergence test: it never
    !!     triggered. Convergence is judged on the age-structure and edge figures
    !!     (scripts/plot_spinup_ageclass_edge.py).

    !! 4.7 Guillaume M. -- MICRO term as a diagnostic. Computed HERE, after the cap, as the
    !!     excess over the prescribed base: that definition holds on BOTH paths, whereas
    !!     edge_src exists only on the regrowth one and is undefined on the legacy tau_rec
    !!     path. On the regrowth path and away from the cap it equals edge_src exactly.
    !!     /!\ Without this field the micro term can only be recovered by subtraction, and
    !!     that subtraction silently breaks as soon as another base term is switched on.
    edge_micro(:) = MAX(zero, edge_length_dyn(:) - edge_length_base(:))

    !! 5. Diagnostics XIOS
    CALL xios_orchidee_send_field("EDGE_LENGTH_DYN", edge_length_dyn)
    ! Guillaume M. -- The scales are output SEPARATELY: without them the decomposition is
    ! not observable, and the departure from the observed edge cannot be attributed term by
    ! term -- which is what the decomposition is for.
    CALL xios_orchidee_send_field("EDGE_MESO",  edge_meso)
    CALL xios_orchidee_send_field("EDGE_MACRO", edge_macro)
    CALL xios_orchidee_send_field("EDGE_HYDRO", edge_hydro)
    CALL xios_orchidee_send_field("EDGE_MICRO", edge_micro)
    ! Guillaume M. -- Management intensity: edge coefficient and patch size actually used
    ! at this step (MODULE_DESIGN_MANAGEMENT_INTENSITY.md).
    CALL xios_orchidee_send_field("ALPHA_HARVEST_EFF", alpha_harvest_use)
    CALL xios_orchidee_send_field("CLEARCUT_AREA",     patch_area_eff)
    ! Guillaume M. -- Target rotation per cell: the ONLY channel by which management
    ! intensity acts on the model dynamics. Output is needed to check that the annual maps
    ! were read, and to compare realised rotation against the target
    ! (scripts/plot_realised_rotation.py).
    IF (ALLOCATED(target_rotation_age)) THEN
       CALL xios_orchidee_send_field("TARGET_ROTATION_AGE", target_rotation_age)
       ! Guillaume M. -- Reference rotation actually used, BEFORE the intensity factor: the
       ! regional map where ROTATION_REF_FILE covers the cell, the namelist scalar elsewhere.
       ! Lets one check that the map was read without going through the harvest cadence.
       rref_out(:,:) = SPREAD(rotation_ref(:), 1, npts)
       IF (ALLOCATED(rotation_ref_map)) THEN
          WHERE (rotation_ref_map(:,:) > zero) rref_out(:,:) = rotation_ref_map(:,:)
       ENDIF
       CALL xios_orchidee_send_field("ROTATION_REF_MAP", rref_out)
       ! Guillaume M. -- Effective terminal-class setpoint seen by the regulators, on every
       ! slot of the group: 1 - a3/R where MAT_A3_ENTRY is active, MAT_TARGET_FRAC(n) elsewhere.
       mtt_out(:,:) = zero
       DO jv = 1, nvm
          IF (.NOT. is_tree(jv)) CYCLE
          IF (nagec_pft(agec_group(jv)) <= 1) CYCLE
          DO ji = 1, npts
             mtt_out(ji,jv) = mat_target_frac_eff(ji, jv, nagec_pft(agec_group(jv)))
          ENDDO
       ENDDO
       CALL xios_orchidee_send_field("MAT_TARGET_TERM", mtt_out)
       ! Guillaume M. -- Online entry ages and the age-based growth rotation; -1 where the
       ! machinery is off or has no sample yet, so the field is always sent.
       mae_out(:,:) = -un
       IF (ALLOCATED(mat_age_entry)) mae_out(:,:) = mat_age_entry(:,:)
       CALL xios_orchidee_send_field("MAT_AGE_ENTRY", mae_out)
       mae_out(:,:) = -un
       IF (ALLOCATED(rotation_growth)) mae_out(:,:) = rotation_growth(:,:)
       CALL xios_orchidee_send_field("ROTATION_GROWTH", mae_out)
    ENDIF
    ! Guillaume M. -- A_open is published as CLASS0_AREA below and nowhere else: once the
    ! landscape term left A_open the former EDGE_REGEN carried the very same values, and two
    ! XIOS fields holding one quantity mislead whoever reads them later. Management-vs-
    ! landscape attribution reads CLASS0_AREA against EDGE_MESO.
    CALL xios_orchidee_send_field("DA_FIRE",    dA_fire)
    CALL xios_orchidee_send_field("DA_STORM",   dA_storm)
    CALL xios_orchidee_send_field("DA_PEST",    dA_pest)
    CALL xios_orchidee_send_field("DA_HARVEST", dA_harvest)

    !! 5b. Guillaume M. -- AED_REGROWTH diagnostics: class-0 stock, graduation flux and the
    !!     management-intensity fraction f_class0 = A_class0 / area_forest. Gated by
    !!     ok_edge_from_age_class, so a flag-off run never references those XIOS fields.
    IF (ok_edge_from_age_class .AND. ALLOCATED(class0_area)) THEN
       ! v2: class0_tot already holds A_open (set in the update step above) -> FRAC_CLASS0 <= 1
       CALL xios_orchidee_send_field("CLASS0_AREA",      class0_tot)
       CALL xios_orchidee_send_field("CLASS0_RECOVERED", class0_recovered)
       IF (ALLOCATED(area_forest)) THEN
          WHERE (area_forest(:) > min_stomate)
             frac_class0_diag(:) = class0_tot(:) / (area_forest(:) * 1.e4_r_std)
          ELSEWHERE
             frac_class0_diag(:) = zero
          ENDWHERE
          CALL xios_orchidee_send_field("FRAC_CLASS0", frac_class0_diag)
       ENDIF
    ENDIF

    !! 5c. Guillaume M. -- MATURITY_CONVEYOR stock class0_veg, summed over groups and bins:
    !!     the grid-cell fraction under regeneration. Gate is n_class0, NOT
    !!     ok_edge_from_age_class -- two distinct conveyors. /!\ While the feed is not
    !!     written the entry bin is zeroed every year in stomate_lpj, so this field is
    !!     IDENTICALLY ZERO by construction; a non-zero value means the feed is connected.
    IF (n_class0 > 0 .AND. ALLOCATED(class0_veg)) THEN
       CALL xios_orchidee_send_field("CLASS0_VEG", &
            SUM(SUM(class0_veg, DIM=3), DIM=2))
    ENDIF

    !! 6. Reset accumulators for next timestep
    dA_fire(:)    = zero
    dA_storm(:)   = zero
    dA_pest(:)    = zero
    dA_harvest(:) = zero

  END SUBROUTINE update_edge_length



!! ================================================================================================
!! FUNCTION   : mat_target_frac_eff
!!
!>\BRIEF        Effective target area share of age class ic within the group of PFT jv, on point ipts.
!!
!! DESCRIPTION: Guillaume M. -- The terminal age class is entered by DIAMETER, so its entry age a3
!!              does not follow the rotation: at equilibrium its share is 1 - a3/R. With
!!              MAT_A3_ENTRY > 0 the terminal setpoint becomes MIN(MAT_A3_FRAC_MAX, MAX(MAT_TARGET_FRAC(n),
!!              1 - a3/target_rotation_age)) and the lower classes are rescaled so the shares sum to
!!              one. Inactive (a3 <= 0, no rotation) = MAT_TARGET_FRAC(ic) untouched, bit-neutral.
!!              See design/MODULE_DESIGN_MAT_TARGET_A3.md.
!! ================================================================================================
  FUNCTION mat_target_frac_eff(ipts, jv, ic)

    IMPLICIT NONE

    INTEGER(i_std), INTENT(in) :: ipts                !! Grid point (local index)
    INTEGER(i_std), INTENT(in) :: jv                  !! Any PFT of the age-class group
    INTEGER(i_std), INTENT(in) :: ic                  !! Age class rank within the group (1..nagec_pft)
    REAL(r_std)                :: mat_target_frac_eff !! Target share of class ic (-)

    INTEGER(i_std) :: jterm, nlast   !! terminal PFT of the group, its class rank
    REAL(r_std)    :: f4, base4      !! effective and namelist terminal setpoints
    REAL(r_std)    :: a3             !! entry age used: online sample or namelist (yr)

    mat_target_frac_eff = mat_target_frac(ic)
    nlast = nagec_pft(agec_group(jv))
    IF (nlast <= 1) RETURN
    jterm = start_index(agec_group(jv)) + nlast - 1
    ! Guillaume M. -- a3 measured by the model itself where MAT_A3_ONLINE has a sample
    ! (age of the area that entered the terminal slot), the namelist parameter otherwise
    ! (seed and fallback). See design/MODULE_DESIGN_MAT_AGE_ENTRY.md.
    a3 = mat_a3_entry(jterm)
    IF (mat_a3_online .AND. ALLOCATED(mat_age_entry)) THEN
       IF (mat_age_entry(ipts,jterm) > zero) a3 = mat_age_entry(ipts,jterm)
    ENDIF
    IF (a3 <= zero) RETURN
    IF (.NOT. ALLOCATED(target_rotation_age)) RETURN
    IF (target_rotation_age(ipts,jterm) <= min_stomate) RETURN

    base4 = mat_target_frac(nlast)
    f4 = MIN(mat_a3_frac_max, MAX(base4, un - a3 / target_rotation_age(ipts,jterm)))
    IF (ic == nlast) THEN
       mat_target_frac_eff = f4
    ELSE IF (base4 < un) THEN
       mat_target_frac_eff = mat_target_frac(ic) * (un - f4) / (un - base4)
    ENDIF

  END FUNCTION mat_target_frac_eff

!! ================================================================================================
!! FUNCTION   : edge_grain_mtc_of
!!
!! DESCRIPTION: Sub-grid fragmentation grain r (m) of a forest PFT, from its MTC.
!!              Used by the meso-scale edge term EL = A*4*f*(1-f)/r. Per MTC rather than
!!              per biome: the MTC already carries the forest type, hence the landscape
!!              structure, and needs no map. Small r = finely divided landscape (broadleaf
!!              in a farmed matrix), large r = coarse blocks (boreal conifer).
!!              CALIBRATED 2026-07-30 against Edge_length_1deg.nc by inverting
!!              r = A*4*f*(1-f)/EL_obs per cell (scripts/calibrate_edge_grain.py, 317 cells
!!              with 0.05<f<0.95). Per-MTC MEDIANS: the mean diverges where f approaches 0
!!              or 1. Acceptance criterion met: 80 % of cells land in 1e2-1e3 m.
!!              /!\ The six MTC agree within 30 % (320-416 m), so the per-MTC structure is
!!              a convention here, not a measured contrast.
!! ================================================================================================
  FUNCTION edge_grain_mtc_of(jv)

    IMPLICIT NONE

    INTEGER(i_std), INTENT(in) :: jv                 !! PFT index
    REAL(r_std)                :: edge_grain_mtc_of  !! grain (m)

    edge_grain_mtc_of = edge_grain_default
    SELECT CASE (pft_to_mtc(jv))
    CASE (4)  ; edge_grain_mtc_of = 416._r_std   ! temperate needleleaf evergreen
    CASE (5)  ; edge_grain_mtc_of = 320._r_std   ! temperate broadleaf evergreen
    CASE (6)  ; edge_grain_mtc_of = 387._r_std   ! temperate broadleaf summergreen
    CASE (7)  ; edge_grain_mtc_of = 414._r_std   ! boreal needleleaf evergreen
    CASE (8)  ; edge_grain_mtc_of = 414._r_std   ! boreal broadleaf summergreen
    CASE (9)  ; edge_grain_mtc_of = 386._r_std   ! boreal needleleaf summergreen (larch)
    ! MTC 2 and 3 (tropical) fall back to edge_grain_default: the EFDA map is pan-EU
    ! only, so they are NOT calibrated. Do not read them as a tropical estimate.
    END SELECT

  END FUNCTION edge_grain_mtc_of


END MODULE stomate_data
