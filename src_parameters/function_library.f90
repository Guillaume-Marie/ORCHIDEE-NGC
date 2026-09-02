! ================================================================================================================================
! MODULE       : function_library
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Collection of functions that are used throughout the ORCHIDEE code
!!
!!\n DESCRIPTION: Collection of modules to : (1) convert one variable into another i.e. basal area
!! to diameter, diamter to tree height, diameter to crown area, etc. (2) ...
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/function_library.f90 $
!! $Date: 2026-02-23 19:54:56 +0100 (lun. 23 févr. 2026) $
!! $Revision: 9393 $
!! \n
!_ ================================================================================================================================

MODULE function_library

  ! modules used:
  USE grid
  USE pft_parameters
  USE dynamic_parameters
  USE constantes_var
  USE mod_orchidee_para_var, ONLY : is_mpi_root, is_omp_root
  USE mod_orchidee_para, ONLY : Set_stdout_file
#ifdef CPP_IEEE_ARITHMETIC
  USE,INTRINSIC :: IEEE_ARITHMETIC, only : IEEE_IS_NAN
#endif
  
  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC calculate_c0_alloc, &
       wood_to_ba_eff, &
       wood_to_ba, &
       cc_kill_to_area, &
       wood_to_height, &
       wood_to_qmheight, &
       wood_to_height_eff, &
       wood_to_tot_volume, &
       wood_to_stand_volume, &
       wood_to_stand_volume_inv, &
       wood_to_dia, &
       wood_to_qmdia,  &
       wood_to_qmdia_up_half, &
       wood_to_circ, &
       wood_to_cn_dia, &
       wood_to_cv, &
       cc_to_biomass, &
       biomass_to_cc, &
       cc_to_lai, &
       biomass_to_lai, &
       lai_to_biomass, &
       Nmax, &
       calculate_rdi, &
       check_vegetation_area, &
       check_area_invariant, &
       check_vegetation_area_point, &
       check_mass_balance, &
       check_mass_balance_point, &
       intermediate_mass_balance_check, &
       check_pixel_area, &
       check_area_change, &
       check_variable_change, &
       check_variable_swap, &
       check_biomass_change, &
       biomass_to_coupled_lai, &
       sort_circ_class_biomass, &
       Arrhenius, &
       Arrhenius_modified, &
       Arrhenius_modified_nodim, &
       get_printlev, &
       ludcmp, &
       swap, &
       imaxloc, &
       mprove, &
       lubksb, &
       polyfit, &
       fm_polynomial_parameter_fitting, &
       calculate_rdi_boundaries, &
       test_param, &
       cc_grid_max_x, &
       weibull_class_dist

  INTERFACE biomass_to_lai
     MODULE PROCEDURE biomass_to_lai_0d, biomass_to_lai_1d
  END INTERFACE

  INTERFACE lai_to_biomass
     MODULE PROCEDURE lai_to_biomass_0d,  lai_to_biomass_1d
  END INTERFACE

  INTERFACE cc_to_biomass
     MODULE PROCEDURE cc_to_biomass_2d, cc_to_biomass_4d
  END INTERFACE cc_to_biomass

  INTERFACE biomass_to_cc
     MODULE PROCEDURE biomass_to_cc_1d, biomass_to_cc_3d
  END INTERFACE biomass_to_cc

  INTERFACE wood_to_stand_volume
     MODULE PROCEDURE wood_to_stand_volume_0d, wood_to_stand_volume_2d
  END INTERFACE

  INTERFACE wood_to_tot_volume
     MODULE PROCEDURE wood_to_tot_volume_0d, wood_to_tot_volume_2d
  END INTERFACE

  INTERFACE cc_to_lai
     MODULE PROCEDURE cc_to_lai_0d, cc_to_lai_1d
  END INTERFACE

  INTERFACE Arrhenius
     MODULE PROCEDURE Arrhenius_0d, Arrhenius_1d
  END INTERFACE

  INTERFACE Arrhenius_modified
     MODULE PROCEDURE Arrhenius_modified_0d, Arrhenius_modified_1d
  END INTERFACE

  CONTAINS


!! ================================================================================================================================
!! FUNCTION    : biomass_to_coupled_lai
!!
!>\BRIEF        Calculate the coupled_LAI based on biomass and veget of each pft
!!
!! DESCRIPTION : Calculates the lai that is coupled with the atmosphere  
!! 
!!             
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::coupled_lai (m2/m2)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION biomass_to_coupled_lai(biomass_leaf, veget, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std)                                       :: pft                    !! PFT number (-)
    REAL(r_std)                                          :: veget                  !! 1-Pgap or the vegetation cover
    REAL(r_std)                                          :: biomass_leaf           !! Biomass of an individual tree within a circ 
                                                                                   !! class @tex $(m^{2} ind^{-1})$ @endtex 

    !! 0.2 Output variables         
    REAL(r_std)                                          :: biomass_to_coupled_lai !! What we get out of this function, depends on the on which
                                                                                   !! leaf biomass is passed to the function;
                                                                                   !! is it at the tree or the stand level. At the tree level
                                                                                   !! it corresponds to the leaf area per tree (m2). At the 
                                                                                   !! stand levels it is the fraction of the LAI that is 
                                                                                   !! assumed to interact  with the atmosphere
                                                                                   !!  @tex $(m^{2} m^{-2})$ @endtex
    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                          :: lai                    !! Leaf area index 
                                                                                   !! @tex $(m^{2} m^{-2})$ @endtex               
!_ ================================================================================================================================
 
    lai = biomass_to_lai(biomass_leaf, pft) 
   
    !+++CHECK+++
    ! In a closed canopy not all the leaves fully interact with the 
    ! atmosphere because leaves can shelter each other. In a more 
    ! open canopy most leaves can interact with the atmosphere.
    ! For the moment we are not clear how to calculate the fraction of
    ! the canopy that is coupled to the atmosphere. Also, based on the
    ! simulated evapotranspiration we need the entire canopy to be
    ! coupled to the atmosphere to obtain simulations which are close
    ! to the observations (Jung's upscaled product).
    IF (veget .GT. min_sechiba) THEN
     
       biomass_to_coupled_lai = lai   

    ELSE
   
       ! There are so many gaps that veget is extremely small. It is fair to
       ! assume that the whole canopy is coupled to the atmosphere.
       biomass_to_coupled_lai = lai
      
    ENDIF

  END FUNCTION biomass_to_coupled_lai
 

!! ================================================================================================================================
!! FUNCTION     : biomass_to_lai_0d
!!
!>\BRIEF        Calculate the LAI based on the leaf biomass
!!
!! DESCRIPTION : Calculates the LAI of a PFT/grid square based on the leaf biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::LAI [m**2 m**{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION biomass_to_lai_0d(leaf_biomass, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std)                                             :: leaf_biomass      !! Biomass of the leaves
                                                                                 !! @tex $(gC m^{-2})$ @endtex 


    !! 0.2 Output variables
                
    REAL(r_std)                                             :: biomass_to_lai_0d !! Leaf area index 
                                                                                 !! @tex $(m^{2} m^{-2})$ @endtex 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                             :: impose_lai         !! LAI read from run.def
!_ ================================================================================================================================
 
    !! 1. Calculate the LAI from the leaf biomass
    IF(sla_dyn) THEN
       biomass_to_lai_0d = log(1.+ ext_coeff_N(pft) * leaf_biomass * slainit(pft))/(ext_coeff_N(pft)) 
    ELSE
       biomass_to_lai_0d = leaf_biomass * sla(pft)
    ENDIF
    
!!$    !+++HACK+++
!!$    ! This is the perfect place to hack the code to make it run with
!!$    ! constant lai
!!$    WRITE(numout,*) 'WARNING: Using fake lai values for testing!'
!!$    biomass_to_lai_1d(:)=3.79052
!!$    !++++++++++

  END FUNCTION biomass_to_lai_0d


!! ================================================================================================================================
!! FUNCTION     : biomass_to_lai_1d
!!
!>\BRIEF        Calculate the LAI based on the leaf biomass
!!
!! DESCRIPTION : Calculates the LAI of a PFT/grid square based on the leaf biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::LAI [m**2 m**{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION biomass_to_lai_1d(leaf_biomass, npts, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std)                                          :: npts
    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(npts)                            :: leaf_biomass      !! Biomass of the leaves
                                                                                 !! @tex $(gC m^{-2})$ @endtex 

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(npts)                            :: biomass_to_lai_1d !! Leaf area index 
                                                                                 !! @tex $(m^{2} m^{-2})$ @endtex 
!_ ================================================================================================================================
 
    !! 1. Calculate the LAI from the leaf biomass
    IF(sla_dyn) THEN
       biomass_to_lai_1d = log(1.+ ext_coeff_N(pft) * leaf_biomass * slainit(pft))/(ext_coeff_N(pft)) 
    ELSE
       biomass_to_lai_1d = leaf_biomass * sla(pft)
    ENDIF

!!$    !+++HACK+++
!!$    ! This is the perfect place to hack the code to make it run with
!!$    ! constant lai
!!$    WRITE(numout,*) 'WARNING: Using fake lai values for testing!'
!!$    biomass_to_lai_1d(:)=3.79052
!!$    !++++++++++

  END FUNCTION biomass_to_lai_1d


!! ================================================================================================================================
!! FUNCTION     : lai_to_biomass_0d
!!
!>\BRIEF        Calculate leaf biomass from lai
!!
!! DESCRIPTION : Calculates leaf biomass (m2 m-2) from leaf biomass (gC m-2)
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::Leaf biomass [gC m^{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION lai_to_biomass_0d(lai, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std)                                             :: lai               !! Leaf Area Index
                                                                                 !! @tex $(m^{2} m^{-2})$ @endtex


    !! 0.2 Output variables
                
    REAL(r_std)                                             :: lai_to_biomass_0d !! Leaf area index
                                                                                 !! @tex $(gC m^{-2})$ @endtex  
                                                                                  
!_ ================================================================================================================================
 
    !! 1. Calculate the LAI from the leaf biomass
    IF(sla_dyn) THEN
       lai_to_biomass_0d = ( 10**(lai*ext_coeff_N(pft)) - 1.) / &
               (ext_coeff_N(pft) * slainit(pft))       
    ELSE
       lai_to_biomass_0d = lai / sla(pft)
    ENDIF

END FUNCTION lai_to_biomass_0d



!! ================================================================================================================================
!! FUNCTION     : lai_to_biomass_1d
!!
!>\BRIEF        Calculate leaf biomass from lai
!!
!! DESCRIPTION : Calculates leaf biomass (m2 m-2) from leaf biomass (gC m-2)
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::Leaf biomass [gC m^{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION lai_to_biomass_1d(lai, npts, pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std)                                          :: npts
    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(npts)                            :: lai               !! Leaf Area Index
                                                                                 !! @tex $(m^{2} m^{-2})$ @endtex


    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(npts)                            :: lai_to_biomass_1d !! Leaf area index
                                                                                 !! @tex $(gC m^{-2})$ @endtex  
                                                                                 
!_ ================================================================================================================================
 
    !! 1. Calculate the LAI from the leaf biomass
    IF(sla_dyn) THEN
       lai_to_biomass_1d = ( 10**(lai*ext_coeff_N(pft)) - 1.) / &
               (ext_coeff_N(pft) * slainit(pft))       
    ELSE
       lai_to_biomass_1d = lai / sla(pft)
    ENDIF

  END FUNCTION lai_to_biomass_1d
  
  
!! ================================================================================================================================
!! FUNCTION     : calculate_c0_alloc
!!
!>\BRIEF        Calculate the baseline root vs sapwood allocation
!!
!! DESCRIPTION : Calculates the baseline root vs sapwood allocation based on the 
!! parameters of the pipe model (hydraulic conductivities) and the
!! turnover of the different components              
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::c0_alloc (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION calculate_c0_alloc(pft, longevity_eff_root, longevity_eff_sap)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                             :: pft                     !! PFT number (-)
    REAL(r_std)                                :: longevity_eff_root      !! Effective longevity for roots (days)
    REAL(r_std)                                :: longevity_eff_sap       !! Effective longevity for sapwood (days)
    
    !! 0.2 Output variables
                
    REAL(r_std)                                :: calculate_c0_alloc  !! quadratic mean height (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                :: sap_density
    REAL(r_std)                                :: qm_dia            !! quadratic mean diameter (m)

!_ ================================================================================================================================
       
    sap_density = deux * pipe_density(pft) / kilo_to_unit
 
    !! 1. Calculate c0_alloc
    IF ( is_tree(pft) ) THEN

       IF (ok_hydrol_arch) THEN
          calculate_c0_alloc = sqrt(k_belowground(pft)/k_sap(pft)*longevity_eff_sap/longevity_eff_root*sapwood_density(pft))
       ELSE
          calculate_c0_alloc = sqrt(k_belowground(pft)/k_sap(pft)*longevity_eff_sap/longevity_eff_root*sap_density)
       ENDIF

    ! Grasses and croplands   
    ELSE

       !+++CHECK+++
       ! Simply copied the same formulation as for trees but note
       ! that the sapwood in trees vs grasses and crops has a very
       ! meaning. In grasses and crops is structural carbon to ensure
       ! that the allocation works. In trees it really is the sapwood
       calculate_c0_alloc = sqrt(k_belowground(pft)/k_sap(pft)*longevity_eff_sap/longevity_eff_root*sap_density)
       !+++++++++++

    ENDIF ! is_tree(j)

 END FUNCTION calculate_c0_alloc


!! ================================================================================================================================
!! FUNCTION     : wood_to_ba_eff
!!
!>\BRIEF        Calculate effective basal area from woody biomass making use of allometric relationships
!!
!! DESCRIPTION :  Calculate basal area of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships. Effective basal area accounts for both above and below ground carbon
!! and is the basis for the application of the rule of Deleuze and Dhote.
!! (i) woodmass = tree_ff * pipe_density*ba*height
!! (ii) height = pipe_tune2 * sqrt(4/pi*ba) ** pipe_tune_3  
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : effective basal area (m2 ind-1)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION wood_to_ba_eff(biomass_temp, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft                  !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp         !! Biomass of an individual tree within a circ 
                                                                                    !! class @tex $(m^{2} ind^{-1})$ @endtex
    REAL(r_std)                                             :: maxheight            !! pipe_tune2 for this pixel and pft

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_ba_eff       !! Effective basal area of an individual tree within a circ
                                                                                    !! class @tex $(m^{2} ind^{-1})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                    !! Index
    REAL(r_std)                                             :: woodmass_ind         !! Woodmass of an individual tree
                                                                                    !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std)                                             :: temp                 !! temporary variable 

!_ ================================================================================================================================
 
 !! 1. Calculate basal area from woodmass

    IF ( is_tree(pft) ) THEN
       
       DO l = 1,ncirc
         
          ! Woodmass of an individual tree
          woodmass_ind = biomass_temp(l,isapabove) + biomass_temp(l,isapbelow) + &
             biomass_temp(l,iheartabove) + biomass_temp(l,iheartbelow)

          temp = woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)
 
          IF (temp.LT.EPSILON(un)) THEN
             ! Expect an overflow when taking an exponent of such a small number
             wood_to_ba_eff(l) = min_stomate
          ELSE
             ! Basal area of that individual (m2 ind-1)  
             wood_to_ba_eff(l) =(pi/4*(woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)) &
                  **(2./pipe_tune3(pft)))**(pipe_tune3(pft)/(pipe_tune3(pft)+2))
             
          ENDIF
                 
       ENDDO

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_ba_eff', &
            'wood_to_ba_eff is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_ba_eff


!! ================================================================================================================================
!! FUNCTION     : wood_to_ba
!!
!>\BRIEF        Calculate basal area from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate basal area of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships given below. Here basal area is defined in line with its classical 
!! forestry meaning.
!! (i) woodmass = tree_ff * pipe_density*ba*height
!! (ii) height = pipe_tune2 * sqrt(4/pi*ba) ** pipe_tune_3  
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : basal area (m2 ind-1)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION wood_to_ba(biomass_temp, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp      !! Biomass of an individual tree within a circ 
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex 
    REAL(r_std)                                             :: maxheight         !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_ba        !! Basal area of an individual tree within a circ
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                 !! Index
    REAL(r_std)                                             :: woodmass_ind      !! Woodmass of an individual tree
                                                                                 !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std)                                             :: temp              !! Temporary variable

!_ ================================================================================================================================
 
 !! 1. Calculate basal area from woodmass

    IF ( is_tree(pft) ) THEN

       DO l = 1,ncirc
             
          ! Woodmass of an individual tree
          woodmass_ind = biomass_temp(l,iheartabove) + biomass_temp(l,isapabove)

          temp = woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)

          IF (temp.LT.EPSILON(un)) THEN
             ! Expect an overflow when taking an exponent of such a small number
             wood_to_ba(l) = min_stomate
          ELSE
             ! Basal area of that individual (m2 ind-1)  
             wood_to_ba(l) = (pi/4*(woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)) & 
                  **(2./pipe_tune3(pft)))**(pipe_tune3(pft)/(pipe_tune3(pft)+2))              
          ENDIF

       ENDDO

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_ba', &
            'wood_to_ba is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_ba


!! ================================================================================================================================
!! SUBROUTINE   : cc_kill_to_area
!!
!>\BRIEF        Convert the number of trees killed per circumference class into a
!!              perturbed forest surface (m2), to harmonise disturbance outputs (windthrow,
!!              pest, fire) with a single basal-area metric. Also used to size the
!!              edge-generating harvest area for the AED feedback (Marie 2026).
!!
!! DESCRIPTION  : circ-class path. For a tree PFT, the killed basal-area fraction
!!               of the stand is computed from the number of trees killed per circ class
!!               (n_kill_cc, supplied directly by the caller), weighted by the
!!               individual basal area (wood_to_ba). That fraction is mapped to a surface
!!               through the prescribed cover fraction and the grid-cell area :
!!                 damaged_area = MIN( Sum(n_kill*ba) / Sum(circ_class_n*ba), 1 )
!!                                * veget_max * cell_area
!!               Capped at veget_max*cell_area : the perturbed area cannot exceed the PFT
!!               forest area on the grid cell.
!!
!! RETURN VALUE : damaged_area (m2), per grid point
!! ================================================================================================================================

  SUBROUTINE cc_kill_to_area(npts, pft, n_kill_cc, circ_class_biomass, circ_class_n, &
                             veget_max, cell_area, maxheight, damaged_area)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)  :: npts                                           !! Domain size (-)
    INTEGER(i_std), INTENT(in)  :: pft                                            !! PFT index (tree) (-)
    REAL(r_std), INTENT(in)     :: n_kill_cc(npts,ncirc)                          !! Trees killed per circ class (ind m-2)
    REAL(r_std), INTENT(in)     :: circ_class_biomass(npts,ncirc,nparts,nelements)!! Biomass of  per circ class(gC ind-1)
    REAL(r_std), INTENT(in)     :: circ_class_n(npts,ncirc)                       !! Stand density per circ class (ind m-2)
    REAL(r_std), INTENT(in)     :: veget_max(npts)                                !! Prescribed PFT cover fraction (-)
    REAL(r_std), INTENT(in)     :: cell_area(npts)                                !! Grid-cell area (m2)
    REAL(r_std), INTENT(in)     :: maxheight(npts)                                !! pipe_tune2 of this pft (m), for wood_to_ba

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)    :: damaged_area(npts)                             !! Perturbed forest surface (m2)

    !! 0.4 Local variables
    INTEGER(i_std)                :: ji, l                                        !! Indices
    REAL(r_std), DIMENSION(ncirc) :: ba                                           !! Basal area per individual tree (m2 ind-1)
    REAL(r_std)                   :: n_kill                                       !! Trees killed in a circ class (ind m-2)
    REAL(r_std)                   :: ba_kill, ba_stand                            !! Killed / stand basal area (m2 m-2)
    REAL(r_std)                   :: frac                                         !! Killed fraction of the stand (-)

!_ ================================================================================================================================

    damaged_area(:) = zero
    IF ( .NOT. is_tree(pft) ) RETURN

    DO ji = 1,npts

       IF ( veget_max(ji) .LE. min_stomate ) CYCLE

       !! Basal area of an individual tree per circ class (m2 ind-1)
       ba(:) = wood_to_ba(circ_class_biomass(ji,:,:,icarbon), pft, maxheight(ji))

       ba_kill  = zero
       ba_stand = zero
       DO l = 1,ncirc
          !! A stand cannot lose more trees than it holds
          n_kill   = MAX(zero, MIN(n_kill_cc(ji,l), circ_class_n(ji,l)))
          ba_kill  = ba_kill  + n_kill             * ba(l)
          ba_stand = ba_stand + circ_class_n(ji,l) * ba(l)
       ENDDO

       IF ( ba_stand .GT. min_stomate ) THEN
          frac = MIN(ba_kill / ba_stand, un)
       ELSE
          frac = zero
       ENDIF

       damaged_area(ji) = frac * veget_max(ji) * cell_area(ji)

    ENDDO

  END SUBROUTINE cc_kill_to_area


!! ================================================================================================================================
!! FUNCTION     : cc_to_biomass
!!
!>\BRIEF        Calculate total biomass from circ_class_biomass
!!
!! DESCRIPTION : circ_class_biomass is expressed in gC tree-1 and thus needs to
!!               be multiplied by the number of trees per circumference class
!!               circ_class_n expressed in tree m-2 to calculate the total
!!               biomass in gC m-2. 
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : total biomass (gC(N) m-2)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION cc_to_biomass_2d(npts,ivm,ccbio,ccind)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: npts              !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                        :: ivm               !! number of the specific PFT (unitless)
    REAL(r_std), DIMENSION(:,:,:)                     :: ccbio             !! circ_class_biomass. Biomass of an individual 
                                                                           !! tree @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:)                         :: ccind             !! circ_class_n. Number of trees per circumference
                                                                           !! class

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(nparts,nelements)          :: cc_to_biomass_2d  !! biomass @tex $(gC(N) m^{2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: iele,icir,ipar    !! Index

!_ ================================================================================================================================
    
    cc_to_biomass_2d(:,:) = zero

    DO iele = 1,nelements

       ! Note that when used in 2d ivm is a single PFT
       IF ( is_tree(ivm) ) THEN

          DO ipar = 1,nparts

                cc_to_biomass_2d(ipar,iele) = &
                     SUM(ccbio(:,ipar,iele) * ccind(:))
                                        
          ENDDO

       ELSE

          ! Grasses and cropland only have one circumference class
          DO ipar = 1,nparts
             
             cc_to_biomass_2d(ipar,iele) = &
                  ccbio(1,ipar,iele) * ccind(1)
             
          END DO
             
       ENDIF

    ENDDO

 END FUNCTION cc_to_biomass_2d

 
!! ================================================================================================================================

 FUNCTION cc_to_biomass_4d(npts,nvm,ccbio,ccind)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: npts               !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                        :: nvm                !! Total number of PFTs (unitless)
    REAL(r_std), DIMENSION(:,:,:,:,:)                 :: ccbio              !! circ_class_biomass. Biomass of an individual 
                                                                            !! tree @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:)                     :: ccind              !! circ_class_n. Number of trees per circumference
                                                                            !! class

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: cc_to_biomass_4d   !! biomass @tex $(gC(N) m^{2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: iele,icir,ipar,ivm !! Index

!_ ================================================================================================================================
    
    cc_to_biomass_4d(:,:,:,:) = zero

    DO iele = 1,nelements

       ! Note that when used in 4d, nvm rather than ivm
       ! should be passed. Now the function will loop over all
       ! PFTs instead of just 1 PFT for the 2d version of this
       ! function. 
       DO ivm = 2,nvm

          IF ( is_tree(ivm) ) THEN

             DO ipar = 1,nparts

                cc_to_biomass_4d(:,ivm,ipar,iele) = &
                     SUM(ccbio(:,ivm,:,ipar,iele) * ccind(:,ivm,:),2)
                
             ENDDO

          ELSE
             
             ! Grasses and cropland only have one circumference class
             DO ipar = 1,nparts

                cc_to_biomass_4d(:,ivm,ipar,iele) = & 
                     ccbio(:,ivm,1,ipar,iele) * ccind(:,ivm,1)

             ENDDO

          ENDIF

       ENDDO

    ENDDO

 END FUNCTION cc_to_biomass_4d


!! ================================================================================================================================
!! FUNCTION     : biomass_to_cc
!!
!>\BRIEF        Calculate circ_class_biomass from biomass
!!
!! DESCRIPTION : The labile and carbres pools are often dealt with at the stand
!!               level for simplicity but then need to be stored at the plant 
!!               level. This is needed because labile and carbres C is lost when
!!               plants are dying. Also, all plant C and N should be stored in a
!!               single variable. 
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : total biomass (gC(N) m-2)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
 FUNCTION biomass_to_cc_1d(bio,ccbio,ccind)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL(r_std)                                       :: bio               !! biomass. Biomass of a specific pool of 
                                                                           !! a whole stand 
                                                                           !! @tex $(gC(N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:)                         :: ccbio             !! circ_class_biomass. Biomass of a specif pool of 
                                                                           !! the individual trees in each circumference class 
                                                                           !! @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:)                         :: ccind             !! circ_class_n. Number of trees per circumference
                                                                           !! class

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                     :: biomass_to_cc_1d  !! biomass @tex $(gC(N) m^{2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: icir              !! Index
    REAL(r_std), DIMENSION(ncirc)                     :: share_ncirc       !! ratio of specific pool
                                                                           !! across circumference classes

!_ ================================================================================================================================
    
    biomass_to_cc_1d(:) = zero

    ! Convert the pools
    IF (SUM(ccbio(:)*ccind(:)) .GT. zero) THEN
       
       ! Proportional to the present pools
       share_ncirc(:) = (ccbio(:)*ccind(:))/SUM(ccbio(:)*ccind(:))

    ELSEIF (SUM(ccind(:)) .GT. zero) THEN

       ! There are no pools at present. Make it proportional
       ! to the number of trees
       share_ncirc(:) = ccind(:)/SUM(ccind(:))

    ELSE 

       ! There is no biomass in circ_class_biomass or individuals
       share_ncirc(:) = zero


    ENDIF

    share_ncirc(1) = zero
    share_ncirc(1) = 1 - SUM(share_ncirc(:))

    ! Distribute bio over the ncirc circumference classes of biomass_to_cc
    ! biomass_to_cc will eventually be assigned to circ_class_biomass
    DO icir = 1,ncirc

       IF (ccind(icir) .GT. zero) THEN

          biomass_to_cc_1d(icir) = bio / ccind(icir) * share_ncirc(icir)

       ELSE

          biomass_to_cc_1d(icir) = zero

       ENDIF

    ENDDO

 END FUNCTION biomass_to_cc_1d

!! ================================================================================================================================
FUNCTION biomass_to_cc_3d(bio,ccbio,ccind,npts,nvm)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts              !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                        :: nvm               !! Total number of PFTs (unitless)
    REAL(r_std),DIMENSION(:,:)                        :: bio               !! biomass. Biomass of a specific pool of 
                                                                           !! a whole stand. dimensions are npts and nvm 
                                                                           !! @tex $(gC(N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:)                     :: ccbio             !! circ_class_biomass. Biomass of a specif pool of 
                                                                           !! the individual trees in each circumference class 
                                                                           !! @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:)                     :: ccind             !! circ_class_n. Number of trees per circumference
                                                                           !! class
   
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: biomass_to_cc_3d  !! biomass @tex $(gC(N) m^{2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: icir,ipts,ivm     !! Index
    REAL(r_std), DIMENSION(ncirc)                     :: share_ncirc       !! ratio of specific pool
                                                                           !! across circumference classes

!_ ================================================================================================================================
    
    biomass_to_cc_3d(:,:,:) = zero

    DO ipts = 1,npts

       DO ivm = 1,nvm

          ! Convert the pools
          IF (is_tree(ivm)) THEN

             ! Forests have ncirc circ classes
             IF (SUM(ccbio(ipts,ivm,:)*ccind(ipts,ivm,:)) .GT. zero) THEN
       
                share_ncirc(:) = (ccbio(ipts,ivm,:)*ccind(ipts,ivm,:)) / &
                     SUM(ccbio(ipts,ivm,:)*ccind(ipts,ivm,:))

             ELSE

                share_ncirc(:) = zero
                
             ENDIF

          ELSE
            
             ! Grass and crops
             share_ncirc(:) = zero

          ENDIF

          ! Make sure that the shares add up to one
          ! and that grass and croplands are properly
          ! dealt with
          share_ncirc(1) = zero
          share_ncirc(1) = 1 - SUM(share_ncirc(:))
                   
          DO icir = 1,ncirc

             IF (ccind(ipts,ivm,icir) .GT. zero) THEN

                biomass_to_cc_3d(ipts,ivm,icir) = bio(ipts,ivm) / &
                     ccind(ipts,ivm,icir) * share_ncirc(icir)

             ELSE

                biomass_to_cc_3d(ipts,ivm,icir) = zero

             ENDIF ! plants are present in this circ_class

          ENDDO ! icir

       ENDDO ! ivm
       
    ENDDO ! ipts

 END FUNCTION biomass_to_cc_3d

!! ================================================================================================================================







!! ================================================================================================================================
!! FUNCTION     : cc_to_lai_0d
!!
!>\BRIEF        Calculate the LAI based on the circ_class_biomass and circ_class_n
!!
!! DESCRIPTION : Calculates the LAI of a PFT/grid square based on the leaf biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::LAI [m**2 m**{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION cc_to_lai_0d(ccbio,ccind,pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: pft                !! Number of the PFT under study (unitless)
    REAL(r_std), DIMENSION(:)                         :: ccbio              !! circ_class_biomass. Biomass of an individual 
                                                                            !! tree @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:)                         :: ccind              !! circ_class_n. Number of trees per circumference
                                                                            !! class

    !! 0.2 Output variables
    REAL(r_std)                                       :: cc_to_lai_0d       !! leaf area index @tex $(m^{2} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: icir               !! Index
    REAL(r_std)                                       :: leaf_biomass       !! Leaf mass for carbon calculated from 
                                                                            !! circ_class_biomass and circ_class_n.
                                                                            !! @tex $(gC m^{-2})$ @endtex
!_ ================================================================================================================================
 
    !! 1. Calculate leaf biomass
    cc_to_lai_0d = zero
    leaf_biomass = zero

    IF ( is_tree(pft) ) THEN

       DO icir = 1,ncirc
                
          leaf_biomass = leaf_biomass + ccbio(icir) * ccind(icir)
                                        
       ENDDO

    ELSE

       ! Grasses and cropland only have one circumference class
       leaf_biomass = ccbio(1) * ccind(1)
        
    ENDIF

    !! 2. Calculate the LAI from the leaf biomass
    IF(sla_dyn) THEN

       cc_to_lai_0d = log(1.+ ext_coeff_N(pft) * leaf_biomass * &
            slainit(pft))/(ext_coeff_N(pft)) 
        
    ELSE
       
       cc_to_lai_0d = leaf_biomass * sla(pft)

    ENDIF

 END FUNCTION cc_to_lai_0d


!! ================================================================================================================================
!! FUNCTION     : cc_to_lai_1d
!!
!>\BRIEF        Calculate the LAI based on the circ_class_biomass and circ_class_n
!!
!! DESCRIPTION : Calculates the LAI of a PFT/grid square based on the leaf biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::LAI [m**2 m**{-2}]
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION cc_to_lai_1d(npts,ccbio,ccind,pft)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts               !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                        :: pft                !! number of PFT under study (unitless)
    REAL(r_std), DIMENSION(:,:,:,:)                   :: ccbio              !! circ_class_biomass. Biomass of an individual 
                                                                            !! tree @tex $(gC(N) ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:)                       :: ccind              !! circ_class_n. Number of trees per circumference
                                                                            !! class

    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(npts)                      :: cc_to_lai_1d       !! leaf area index @tex $(m^{2} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: icir               !! Index
    REAL(r_std), DIMENSION(npts)                      :: leaf_biomass       !! Leaf mass for carbon calculated from 
                                                                            !! circ_class_biomass and circ_class_n.
                                                                            !! @tex $(gC m^{-2})$ @endtex
!_ ================================================================================================================================
 
    !! 1. Calculate leaf biomass for a specific PFT
    cc_to_lai_1d(:) = zero
    leaf_biomass(:) = zero

    ! Note that when used in 2d ivm is a single PFT
    IF ( is_tree(pft) ) THEN

       DO icir = 1,ncirc
                
          leaf_biomass(:) = leaf_biomass(:) + &
               ccbio(:,icir,ileaf,icarbon) * ccind(:,icir)
                                        
       ENDDO

    ELSE

       ! Grasses and cropland only have one circumference class
             
        leaf_biomass(:) = ccbio(:,1,ileaf,icarbon) * ccind(:,1)

     ENDIF

     !! 2. Calculate the LAI from the leaf biomass
     IF(sla_dyn) THEN

       cc_to_lai_1d(:) = log(1.+ ext_coeff_N(pft)* leaf_biomass(:) * &
            slainit(pft))/(ext_coeff_N(pft)) 

    ELSE

       cc_to_lai_1d(:) = leaf_biomass(:) * sla(pft)

    ENDIF

 END FUNCTION cc_to_lai_1d


!! ================================================================================================================================
!! FUNCTION     : wood_to_qmheight
!!
!>\BRIEF        Calculate the quadratic mean height from the biomass
!!
!! DESCRIPTION : Calculates the quadratic mean height from the biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::qm_height (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_qmheight(biomass_temp, ind, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                             :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(ncirc,nparts)       :: biomass_temp      !! Biomass of the leaves @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)              :: ind               !! Number of individuals @tex $(m^{-2})$ @endtex
    REAL(r_std)                                :: maxheight         !! pipe_tune2 for this pixel and pft

    !! 0.2 Output variables
                
    REAL(r_std)                                :: wood_to_qmheight  !! quadratic mean height (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(ncirc)              :: circ_class_ba     !! basal area for each circ_class @tex $(m^{2})$ @endtex
    REAL(r_std)                                :: qm_dia            !! quadratic mean diameter (m)

!_ ================================================================================================================================
 
    !! 1. Calculate qm_height from the biomass
    IF ( is_tree(pft) ) THEN

       ! Basal area at the tree level (m2 tree-1)
       circ_class_ba(:) = wood_to_ba(biomass_temp(:,:),pft,maxheight) 
       
       IF (SUM(ind(:)) .NE. zero) THEN

          qm_dia = SQRT( 4/pi*SUM(circ_class_ba(:)*ind(:))/SUM(ind(:)) )
          
       ELSE
          
          qm_dia = zero

       ENDIF
       
       wood_to_qmheight = maxheight*(qm_dia**pipe_tune3(pft))
       
    ELSE

       ! Grasses and croplands
       ! Calculate height as a function of the leaf biomass. Ensure that the 
       ! grasslands have a roughness length during the winter. Therefore use a
       ! lower threshold.
       IF (SUM(ind(:)) .NE. zero) THEN
          wood_to_qmheight = MAX(0.2, biomass_to_lai(biomass_temp(1,ileaf),pft)&
                * ind(1) * lai_to_height(pft)) 

          !wood_to_qmheight = MAX(0.2,biomass_temp(1,ileaf)*ind(1)*sla(pft)*lai_to_height(pft))
       ELSE

          wood_to_qmheight = zero

       ENDIF

    ENDIF ! is_tree(j)

 END FUNCTION wood_to_qmheight


!! ================================================================================================================================
!! FUNCTION     : wood_to_qmdia
!!
!>\BRIEF        Calculate the quadratic mean diameter from the biomass
!!
!! DESCRIPTION : Calculates the quadratic mean diameter from the aboveground biomass
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::qm_dia (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_qmdia(biomass_temp, ind, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                             :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(ncirc,nparts)       :: biomass_temp      !! Biomass of the leaves @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)              :: ind               !! Number of individuals @tex $(m^{-2})$ @endtex
    REAL(r_std)                                :: maxheight         !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std)                                :: wood_to_qmdia     !! quadratic mean diameter (m) 

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(ncirc)              :: circ_class_ba     !! basal area for each circ_class @tex $(m^{2})$ @endtex

!_ ================================================================================================================================
 
    !! 1. Calculate qm_dia from the biomass
    IF ( is_tree(pft) ) THEN

       ! Basal area at the tree level (m2 tree-1)
       circ_class_ba(:) = wood_to_ba(biomass_temp(:,:),pft,maxheight) 
       
       IF (SUM(ind(:)) .NE. zero) THEN

          wood_to_qmdia = SQRT( 4/pi*SUM(circ_class_ba(:)*ind(:))/SUM(ind(:)) )
          
       ELSE
          
          wood_to_qmdia = zero

       ENDIF


    ! Grasses and croplands   
    ELSE

       wood_to_qmdia = zero   
       
    ENDIF ! is_tree(pft)

 END FUNCTION wood_to_qmdia

!! ================================================================================================================================
!! FUNCTION     : wood_to_qmdia_up_half
!!
!>\BRIEF        Calculate the quadratic mean diameter from the biomass of the upper half of the diameter classes
!!
!! DESCRIPTION : -
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::qm_diad_up_half (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

 FUNCTION wood_to_qmdia_up_half(biomass_temp, ind, pft, maxheight)

   !! 0. Variable and parameter declaration
  
      !! 0.1 Input variables
  
      INTEGER(i_std)                             :: pft                    !! PFT number (-)
      REAL(r_std), DIMENSION(ncirc,nparts)       :: biomass_temp           !! Biomass of the leaves @tex $(gC m^{-2})$ @endtex 
      REAL(r_std), DIMENSION(ncirc)              :: ind                    !! Number of individuals @tex $(m^{-2})$ @endtex
      REAL(r_std)                                :: maxheight              !! pipe_tune2 for this pixel and pft
      
      !! 0.2 Output variables
                  
      REAL(r_std)                                :: wood_to_qmdia_up_half  !! quadratic mean diameter (m) 

  
      !! 0.3 Modified variables
  
      !! 0.4 Local variables
      INTEGER(i_std)                             :: med                    !! circumference class at which dominante diameter is calculated (-)
      REAL(r_std), DIMENSION(ncirc)              :: circ_class_ba          !! basal area for each circ_class @tex $(m^{2})$ @endtex
  
  !_ ================================================================================================================================
   
      !! 1. Calculate qm_dia from the biomass
      IF ( is_tree(pft) ) THEN
  
         ! Basal area at the tree level (m2 tree-1)
         circ_class_ba(:) = wood_to_ba(biomass_temp(:,:),pft,maxheight) 

         med = ceiling(ncirc/deux)

         IF (SUM(ind(med:ncirc)) .NE. zero) THEN
  
            wood_to_qmdia_up_half = SQRT( 4/pi*SUM(circ_class_ba(med:ncirc)*ind(med:ncirc))/SUM(ind(med:ncirc)))
            
         ELSE
            
            wood_to_qmdia_up_half = zero
  
         ENDIF
  
  
      ! Grasses and croplands   
      ELSE
  
         wood_to_qmdia_up_half = zero   
         
      ENDIF ! is_tree(pft)
  
   END FUNCTION wood_to_qmdia_up_half


!! ================================================================================================================================
!! FUNCTION     : wood_to_stand_volume_0d
!!
!>\BRIEF        This allometric function computes stand volume (aboveground) as a function of 
!! biomass at stand scale. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : wood_to_stand_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_stand_volume_0d(npts,ccbio,ccind,pft,branch_ratio,inc_branches)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                :: npts           !! Domain size (unitless)
    REAL(r_std), INTENT(in), DIMENSION(:,:,:) :: ccbio          !! circ_class_biomass or biomass of an 
                                                                !! individual tree in each circumference 
                                                                !! class @tex $(gC tree^{-1})$ @endtex  
    REAL(r_std), INTENT(in), DIMENSION(:)     :: ccind          !! Number of trees per circumference class
    REAL(r_std), INTENT(in)                   :: branch_ratio   !! Branch ratio of sap and heartwood biomass
                                                                !! unitless
    INTEGER(i_std), INTENT(in)                :: pft            !! Plant functional type (unitless)
    INTEGER(i_std), INTENT(in)                :: inc_branches   !! Include the branches in the volume calculation?
                                                                !! 0: exclude the branches from the volume calculation 
                                                                !! (thus correct the biomass for the branch ratio)
                                                                !! 1: include the branches in the volume calculation  
                                                                !! (thus use all aboveground biomass)
                                    
    

    !! 0.2 Output variables

    REAL(r_std)                               :: wood_to_stand_volume_0d !! The volume of aboveground wood per square meter 
                                                                !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nparts,nelements)  :: biomass        !! Total biomass at the stand level
                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                               :: woody_biomass  !! Woody biomass at the stand level 
                                                                !! @tex $(gC m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. Volume to biomass

    ! Calculate stand biomass from circ_class_biomass. In cc_to_biomass_2d
    ! which is used here, npts is not used but we need to pass it to have
    ! a consistent call structure for cc_to_biomass_2d and cc_to_biomass_4d
    biomass = cc_to_biomass(npts,pft,ccbio,ccind)

    ! Woody biomass used in the calculation
    IF (inc_branches .EQ. 0) THEN

       ! The branches are excluded from the calculated volume
       woody_biomass = (biomass(isapabove,icarbon)+biomass(iheartabove,icarbon)) * &
            (un - branch_ratio)

    ELSEIF (inc_branches .EQ. 1) THEN

       ! The branches are included in the calculated volume
       woody_biomass = (biomass(isapabove,icarbon)+biomass(iheartabove,icarbon))

    ELSE

       CALL ipslerr_p(3,'ERROR in wood_to_volume',&
            'Specify how the branches should be accounted for','','')

    ENDIF

    ! Wood volume expressed in m**3 / m**2
    wood_to_stand_volume_0d = woody_biomass/(pipe_density(pft))

  END FUNCTION wood_to_stand_volume_0d


!! ================================================================================================================================
!! FUNCTION     : wood_to_stand_volume_2d
!!
!>\BRIEF        This allometric function computes stand volume (aboveground) as a function of 
!! biomass at stand scale. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : wood_to_stand_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_stand_volume_2d(npts,ccbio,ccind,branch_ratio,inc_branches)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts           !! Domain size (unitless)
    REAL(r_std), INTENT(in), DIMENSION(:,:,:,:,:)     :: ccbio          !! circ_class_biomass or biomass of an 
                                                                        !! individual tree in each circumference 
                                                                        !! class @tex $(gC tree^{-1})$ @endtex  
    REAL(r_std), INTENT(in), DIMENSION(:,:,:)         :: ccind          !! Number of trees per circumference class
    REAL(r_std), INTENT(in), DIMENSION(:)             :: branch_ratio   !! Branch ratio of sap and heartwood biomass
                                                                        !! unitless
  
    INTEGER(i_std), INTENT(in)                        :: inc_branches   !! Include the branches in the volume calculation?
                                                                        !! 0: exclude the branches from the volume calculation 
                                                                        !! (thus correct the biomass for the branch ratio)
                                                                        !! 1: include the branches in the volume calculation  
                                                                        !! (thus use all aboveground biomass)
                                    
    

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm)                  :: wood_to_stand_volume_2d !! The volume of wood per square meter 
                                                                        !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: ipts, ivm      !! indices     
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: biomass        !! Total biomass at the stand level
                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: woody_biomass  !! Woody biomass at the stand level 
                                                                        !! @tex $(gC m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. biomass to volume
    wood_to_stand_volume_2d(:,:) = zero

    ! Calculate stand biomass from circ_class_biomass.
    biomass(:,:,:,:) = cc_to_biomass(npts,nvm,ccbio,ccind)

    ! Woody biomass used in the calculation
    IF (inc_branches .EQ. 0) THEN

       ! The branches are excluded from the calculated volume
       DO ivm = 2,nvm
          woody_biomass(:,ivm) = (biomass(:,ivm,isapabove,icarbon) + &
               biomass(:,ivm,iheartabove,icarbon)) * &
               (un - branch_ratio(ivm))
       END DO

    ELSEIF (inc_branches .EQ. 1) THEN

       ! The branches are included in the calculated volume
       woody_biomass(:,:) = (biomass(:,:,isapabove,icarbon) + &
            biomass(:,:,iheartabove,icarbon))

    ELSE

       CALL ipslerr_p(3,'ERROR in wood_to_volume', &
            'Specify how the branches should be accounted for','','')

    ENDIF

    ! Only calculate for forests
    DO ivm = 2,nvm

       IF ( is_tree(ivm) ) THEN

          ! Wood volume expressed in m**3 / m**2
          wood_to_stand_volume_2d(:,ivm) = woody_biomass(:,ivm)/(pipe_density(ivm))

       END IF

    END DO

  END FUNCTION wood_to_stand_volume_2d


!! ================================================================================================================================
!! FUNCTION     : wood_to_tot_volume_0d
!!
!>\BRIEF        This allometric function computes total volume (above and belowground) as a function of 
!! biomass at stand scale. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : wood_to_tot_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_tot_volume_0d(npts,ccbio,ccind,pft,branch_ratio,inc_branches)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                :: npts           !! Domain size (unitless)
    REAL(r_std), INTENT(in), DIMENSION(:,:,:) :: ccbio          !! circ_class_biomass or biomass of an 
                                                                !! individual tree in each circumference 
                                                                !! class @tex $(gC tree^{-1})$ @endtex  
    REAL(r_std), INTENT(in), DIMENSION(:)     :: ccind          !! Number of trees per circumference class
    REAL(r_std), INTENT(in)                   :: branch_ratio   !! Branch ratio of sap and heartwood biomass
                                                                !! unitless
    INTEGER(i_std), INTENT(in)                :: pft            !! Plant functional type (unitless)
    INTEGER(i_std), INTENT(in)                :: inc_branches   !! Include the branches in the volume calculation?
                                                                !! 0: exclude the branches from the volume calculation 
                                                                !! (thus correct the biomass for the branch ratio)
                                                                !! 1: include the branches in the volume calculation  
                                                                !! (thus use all aboveground biomass)
                                    
    

    !! 0.2 Output variables

    REAL(r_std)                               :: wood_to_tot_volume_0d !! The total volume (above and belowground) of wood per square meter 
                                                                !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(nparts,nelements)  :: biomass        !! Total biomass at the stand level
                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                               :: woody_biomass  !! Woody biomass at the stand level 
                                                                !! @tex $(gC m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. Volume to biomass

    ! Calculate stand biomass from circ_class_biomass. In cc_to_biomass_2d
    ! which is used here, npts is not used but we need to pass it to have
    ! a consistent call structure for cc_to_biomass_2d and cc_to_biomass_4d
    biomass = cc_to_biomass(npts,pft,ccbio,ccind)

    ! Woody biomass used in the calculation
    IF (inc_branches .EQ. 0) THEN

       ! The branches are excluded from the calculated volume
       woody_biomass = (biomass(isapabove,icarbon) + biomass(iheartabove,icarbon) + &
            biomass(isapbelow,icarbon) + biomass(iheartbelow,icarbon)) * &
            (un - branch_ratio)

    ELSEIF (inc_branches .EQ. 1) THEN

       ! The branches are included in the calculated volume
       woody_biomass = biomass(isapabove,icarbon) + biomass(iheartabove,icarbon) + &
            biomass(isapbelow,icarbon) + biomass(iheartbelow,icarbon)

    ELSE

       CALL ipslerr_p(3,'ERROR in wood_to_volume',&
            'Specify how the branches should be accounted for','','')

    ENDIF

    ! Wood volume expressed in m**3 / m**2
    wood_to_tot_volume_0d = woody_biomass/(pipe_density(pft))

  END FUNCTION wood_to_tot_volume_0d


!! ================================================================================================================================
!! FUNCTION     : wood_to_tot_volume_2d
!!
!>\BRIEF        This allometric function computes total volume (above and belowground) as a function of 
!! biomass at stand scale. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : wood_to_tot_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_tot_volume_2d(npts,ccbio,ccind,branch_ratio,inc_branches)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts           !! Domain size (unitless)
    REAL(r_std), INTENT(in), DIMENSION(:,:,:,:,:)     :: ccbio          !! circ_class_biomass or biomass of an 
                                                                        !! individual tree in each circumference 
                                                                        !! class @tex $(gC tree^{-1})$ @endtex  
    REAL(r_std), INTENT(in), DIMENSION(:,:,:)         :: ccind          !! Number of trees per circumference class
    REAL(r_std), INTENT(in), DIMENSION(:)             :: branch_ratio   !! Branch ratio of sap and heartwood biomass
                                                                        !! unitless
  
    INTEGER(i_std), INTENT(in)                        :: inc_branches   !! Include the branches in the volume calculation?
                                                                        !! 0: exclude the branches from the volume calculation 
                                                                        !! (thus correct the biomass for the branch ratio)
                                                                        !! 1: include the branches in the volume calculation  
                                                                        !! (thus use all aboveground biomass)
                                    
    

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm)                  :: wood_to_tot_volume_2d !! The total volume (above and belowground) of wood per square meter 
                                                                        !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: ipts, ivm      !! indices     
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: biomass        !! Total biomass at the stand level
                                                                        !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: woody_biomass  !! Woody biomass at the stand level 
                                                                        !! @tex $(gC m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. biomass to volume
    wood_to_tot_volume_2d(:,:) = zero

    ! Calculate stand biomass from circ_class_biomass.
    biomass(:,:,:,:) = cc_to_biomass(npts,nvm,ccbio,ccind)

    ! Woody biomass used in the calculation
    IF (inc_branches .EQ. 0) THEN

       ! The branches are excluded from the calculated volume
       DO ivm = 2,nvm
          woody_biomass(:,ivm) = (biomass(:,ivm,isapabove,icarbon) + &
               biomass(:,ivm,iheartabove,icarbon) + &
               biomass(:,ivm,isapbelow,icarbon) + &
               biomass(:,ivm,iheartbelow,icarbon)) * &
               (un - branch_ratio(ivm))
       END DO

    ELSEIF (inc_branches .EQ. 1) THEN

       ! The branches are included in the calculated volume
       woody_biomass(:,:) = biomass(:,:,isapabove,icarbon) + &
            biomass(:,:,iheartabove,icarbon) + &
            biomass(:,:,isapbelow,icarbon) + &
            biomass(:,:,iheartbelow,icarbon)

    ELSE

       CALL ipslerr_p(3,'ERROR in wood_to_volume', &
            'Specify how the branches should be accounted for','','')

    ENDIF

    ! Only calculate for forests
    DO ivm = 2,nvm

       IF ( is_tree(ivm) ) THEN

          ! Wood volume expressed in m**3 / m**2
          wood_to_tot_volume_2d(:,ivm) = woody_biomass(:,ivm)/(pipe_density(ivm))

       END IF

    END DO

  END FUNCTION wood_to_tot_volume_2d


!! ================================================================================================================================
!! FUNCTION     : wood_to_stand_volume_inv
!!
!>\BRIEF        This allometric function computes volume as a function of 
!! biomass at stand scale. It only takes trees into account for the diameter exceeds
!! a user defined threshold (e.g. 7cm). It thus better matches what is being done in 
!! forest inventories. Volume \f$(m^3 m^{-2}) = f(biomass (gC m^{-2}))\f$
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!! 
!! RETURN VALUE : biomass_to_volume
!!
!! REFERENCE(S)	: See above, module description.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  FUNCTION wood_to_stand_volume_inv(npts,ccbio,ccind,branch_ratio,inc_branches,maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts               !! Domain size (unitless)
    REAL(r_std), INTENT(in), DIMENSION(:,:,:,:,:)     :: ccbio              !! circ_class_biomass or biomass of an 
                                                                            !! individual tree in each circumference 
                                                                            !! class @tex $(gC tree^{-1})$ @endtex  
    REAL(r_std), INTENT(in), DIMENSION(:,:,:)         :: ccind              !! Number of trees per circumference class
    REAL(r_std), INTENT(in), DIMENSION(:)             :: branch_ratio       !! Branch ratio of sap and heartwood biomass
                                                                            !! unitless
  
    INTEGER(i_std), INTENT(in)                        :: inc_branches       !! Include the branches in the volume calculation?
                                                                            !! 0: exclude the branches from the volume calculation 
                                                                            !! (thus correct the biomass for the branch ratio)
                                                                            !! 1: include the branches in the volume calculation  
                                                                            !! (thus use all aboveground biomass)
    REAL(r_std)                                       :: maxheight          !! pipe_tune2 for this pixel and pft             
    

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm)                  :: wood_to_stand_volume_inv !! The volume of wood per square meter 
                                                                            !!  @tex $(m^3 m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: ipts, ivm          !! indices     
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: biomass            !! Total biomass at the stand level
                                                                            !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: woody_biomass      !! Woody biomass at the stand level 
                                                                            !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: circ_n             !! temporary circ_class_n for estimates to 
                                                                            !! compare against inventory data
    REAL(r_std), DIMENSION(ncirc)                     :: circ_dia           !! diameter for each circumference class (m)

!_ ================================================================================================================================

    ! Initialize
    wood_to_stand_volume_inv(:,:) = zero

    ! Calculate wood volume as it is done in forest inventories
    DO ivm = 2,nvm

       IF ( is_tree(ivm) ) THEN
          
          DO ipts = 1,npts
             
             circ_n(:) = ccind(ipts,ivm,:)
             circ_dia(:) = wood_to_dia(ccbio(ipts,ivm,:,:,icarbon),ivm,maxheight)

             ! Mask diameter classes that do need meet the
             ! inventory diameter criterion
             WHERE (circ_dia(:) .lt. dia_thresh_inv(ivm))
                circ_n(:) = zero
             ENDWHERE

             ! Calculate wood volume
             IF (SUM(circ_n(:)).EQ.zero) THEN
                wood_to_stand_volume_inv(ipts,ivm) = zero
             ELSE
                ! Calculate the wood volume but only account for
                ! diameters > threshold.
                wood_to_stand_volume_inv(ipts,ivm) = wood_to_stand_volume(npts, &
                     ccbio(ipts,ivm,:,:,:), circ_n(:),ivm, &
                     branch_ratio(ivm),0)
             END IF

          ENDDO ! ipts

       ENDIF ! is_tree

    ENDDO ! ivm

  END FUNCTION wood_to_stand_volume_inv


!! ================================================================================================================================
!! FUNCTION     : wood_to_height
!!
!>\BRIEF        Calculate tree height from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate height of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships. This is the height used in forestry and for calculating the aerodynamic
!! interactions 
!! (i) height(:) = pipe_tune2(j)*(4/pi*ba(:))**(pipe_tune3(j)/2) and
!! (ii) woodmass_ind = tree_ff*pipe_density*ba*height   
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : height (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_height(biomass_temp, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft                !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp       !! Biomass of an individual tree within a circ 
                                                                                  !! class @tex $(m^{2} ind^{-1})$ @endtex 
    REAL(r_std)                                             :: maxheight          !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_height     !! Height of an individual tree within a circ
                                                                                  !! class (m)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                  !! Index
    REAL(r_std)                                             :: woodmass_ind       !! Woodmass of an individual tree
                                                              !! @tex $(gC ind{-1})$ @endtex
    REAL(r_std)                                             :: temp               !! Temporary variable
!_ ================================================================================================================================
 
 !! 1. Calculate height from woodmass

    IF ( is_tree(pft) ) THEN
       
       DO l = 1,ncirc
         
          ! Woodmass of an individual tree (only the aboveground component)
          woodmass_ind = biomass_temp(l,isapabove) + biomass_temp(l,iheartabove)

          temp = (woodmass_ind/(tree_ff(pft)*pipe_density(pft))*4/pi)

          IF (temp.LT.EPSILON(un)) THEN
             ! Expect an overflow when taking an exponent of such a small number
             wood_to_height(l) = min_stomate
          ELSE
             ! Height of that individual
             wood_to_height(l) = (((woodmass_ind/(tree_ff(pft)*&
                  pipe_density(pft))*4/pi)**(pipe_tune3(pft)/2))*&
                  maxheight)**(1/(pipe_tune3(pft)/2+1))
          ENDIF
            
       ENDDO

    ELSE

       WRITE(numout,*) 'ERROR: wood_to_height, pft ',pft
       CALL ipslerr_p (3,'wood_to_height', &
            'wood_to_height is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_height


!! ================================================================================================================================
!! FUNCTION     : wood_to_height_eff
!!
!>\BRIEF        Calculate the effective tree height from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate the effective height of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships. Effective height makes use of both above and belowground biomass and is
!! used in the calculation of the allocation according to deleuze and dhote, hydraulic architecture because also
!! the height of the belowground part should be included.
!! (i) height(:) = pipe_tune2(j)*(4/pi*ba(:))**(pipe_tune3(j)/2) and
!! (ii) woodmass_ind = tree_ff*pipe_density*ba*height   
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : height (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_height_eff(biomass_temp, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft                !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp       !! Biomass of an individual tree within a circ 
                                                                                  !! class @tex $(m^{2} ind^{-1})$ @endtex 
    REAL(r_std)                                             :: maxheight          !! pipe_tune2 for this pixel and pft

    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_height_eff !! Effective height of an individual tree within a circ
                                                                                  !! class (m)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                  !! Index
    REAL(r_std)                                             :: woodmass_ind       !! Woodmass of an individual tree
                                                                                  !! @tex $(gC ind{-1})$ @endtex
    REAL(r_std)                                             :: temp               !! Temporary variable                                                                              
!_ ================================================================================================================================
 
 !! 1. Calculate height from woodmass

    IF ( is_tree(pft) ) THEN
       
       DO l = 1,ncirc
         
          ! Woodmass of an individual tree. Both above and belowground biomass are used
          woodmass_ind = biomass_temp(l,isapabove) + biomass_temp(l,isapbelow) + &
             biomass_temp(l,iheartabove) + biomass_temp(l,iheartbelow)

          temp = (woodmass_ind/(tree_ff(pft)*pipe_density(pft))*4/pi)
          
          IF (temp.LT.EPSILON(un)) THEN
             ! Expect an overflow when taking an exponent of such a small number
             wood_to_height_eff(l) = min_stomate
          ELSE
             ! Height of that individual
             wood_to_height_eff(l) = (((woodmass_ind/(tree_ff(pft)*pipe_density(pft))*4/pi)**(pipe_tune3(pft)/2))&
                  *maxheight)**(1/(pipe_tune3(pft)/2+1))
          ENDIF
          
       ENDDO

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_height_eff', &
            'wood_to_height_eff is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_height_eff


!! ================================================================================================================================
!! FUNCTION     : wood_to_dia
!!
!>\BRIEF        Calculate diameter from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate diameter of an individual tree from the woody biomass of that tree making 
!! use of allometric relationships. Makes only use of the aboveground biomass and relates to the
!! typical forestry diameter (but not normalized at 1.3 m)
!! (i) woodmass_ind = tree_ff * pipe_density * height * pi/4*dia**2
!! (ii) height = pipe_tune2 * dia * pipe_tune3 
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : diameter (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_dia(biomass_temp, pft, maxheight)


 !! 0. Variable and parameter declaration


    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp      !! Biomass of an individual tree within a circ 
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex 
    REAL(r_std)                                             :: maxheight         !! pipe_tune2 for this pixel and pft

    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_dia       !! Diameter of an individual tree within a circ
                                                                                 !! class (m)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                 !! Index
    REAL(r_std)                                             :: woodmass_ind      !! Woodmass of an individual tree
                                                                                 !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std)                                             :: temp              !! temporary variable
!_ ================================================================================================================================
 
 !! 1. Calculate basal area from woodmass

    IF ( is_tree(pft) ) THEN
       
       DO l = 1,ncirc
         
          ! Woodmass of an individual tree
          woodmass_ind = biomass_temp(l,isapabove) + biomass_temp(l,iheartabove)

          temp = 4/pi*woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)

          IF (temp.LT.EPSILON(un)) THEN
             ! Expect an overflow when taking an exponent of such a small number
             wood_to_dia(l) = min_stomate
          ELSE
             ! Basal area of that individual (m2 ind-1)  
             wood_to_dia(l) = (4/pi*woodmass_ind/(tree_ff(pft)*pipe_density(pft)*maxheight)) ** &
                  (1./(2+pipe_tune3(pft)))
          ENDIF   
       ENDDO

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_dia', &
            'wood_to_dia is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

 END FUNCTION wood_to_dia


!! ================================================================================================================================
!! FUNCTION     : wood_to_circ
!!
!>\BRIEF        Calculate circumference from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : All this does it computer the diameter using a different routine, and then
!!               convert that into a circumference.  
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : circumference (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_circ(biomass_temp, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp      !! Biomass of an individual tree within a circ 
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex 
    REAL(r_std)                                             :: maxheight         !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_circ      !! Circumference of an individual tree within a circ
                                                                                 !! class (m)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================
 
    !! 1. Calculate diameter from woodmass

    wood_to_circ(:)=val_exp

    wood_to_circ(:)=wood_to_dia(biomass_temp(:,:),pft,maxheight)

    ! convert to a circumference (m)
    wood_to_circ(:) =  wood_to_circ(:)*pi

 END FUNCTION wood_to_circ


!! ================================================================================================================================
! SUBROUTINE    : wood_to_cn_dia
!!
!>\BRIEF        Calculate horizontal crown diameter from woody biomass making use of allometric relationships
!!
!! DESCRIPTION : Calculate tree height and apply the reduction coefficient crown_to_height and crown_vertohor_dia 
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : horizontal crown diameter (m)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 SUBROUTINE wood_to_cn_dia(biomass_temp, ind, pft, dia_hor, dia_ver, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft                !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp       !! Biomass of an individual tree within a circ 
                                                                                  !! class @tex $(m^{2} ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:)                               :: ind                !! Number of individuals within a circ class
                                                                                  !! @tex $(ind m^{2})$ @endtex
    REAL(r_std)                                             :: maxheight          !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: dia_hor            !! Horizontal crown diameter 
                                                                                  !! for each individual tree within a circ
                                                                                  !! class (m)
    REAL(r_std), DIMENSION(ncirc)                           :: dia_ver             !! Vertical crown diameter (in othe words, crown 
                                                                                  !! depth) for each individual tree within a circ
                                                                                  !! class (m)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: icir               !! Index
    REAL(r_std)                                             :: woodmass_ind       !! Woodmass of an individual tree
                                                                                  !! @tex $(gC ind{-1})$ @endtex
    REAL(r_std)                                             :: temp               !! Temporary variable                                                                              
    REAL(r_std)                                             :: cn_vol             !! Crown volume of a individual tree within a circ class (m)
    REAL(r_std)                                             :: total_cn_vol       !! Total crown volume of the stand across all circ classes (m)
    REAL(r_std)                                             :: max_height         !! Volume of the entire canopy space (surface area * tree height 
                                                                                  !! but normalized for surface area) (m)
    REAL(r_std)                                             :: adjusted_crown_vertohor_dia !! calculated ratio between vertical and horizontal 
                                                                                  !! crown diameter. The calculated value guarantees that the 
                                                                                  !! canopys space is filled for 74% (given as ::crown_packing)
    REAL(r_std)                                             :: adjusted_crown_to_height !! Calculate ratio between tree height and vertical
                                                                                  !! crown diameter. The calculated value guarantees that the 
                                                                                  !! canopys space is filled for 74% (given as ::crown_packing)

!_ ================================================================================================================================
          
 !! 1. Calculate horizontal crown diameter from woodmass

    IF ( is_tree(pft) ) THEN
       
       ! Initialize
       total_cn_vol=zero

       ! Use wood mass to calculate the prescribed vertical and horizontal crown 
       ! diameters
       DO icir = 1,ncirc
         
          ! Woodmass of an individual tree (only the aboveground component)
          woodmass_ind = biomass_temp(icir,isapabove) + biomass_temp(icir,iheartabove)

          temp = (woodmass_ind/(tree_ff(pft)*pipe_density(pft))*4/pi)

          IF (temp.LT.EPSILON(un)) THEN

             ! Expect an overflow when taking an exponent of such a small number
             dia_hor(icir) = min_stomate
             dia_ver(icir) = min_stomate

          ELSE
             
             ! Horizonatal crown diameter of an individual based on the height
             ! and prescribed ratios.  
             dia_hor(icir) = crown_vertohor_dia(pft) * crown_to_height(pft) * &
                  (((woodmass_ind/(tree_ff(pft)*&
                  pipe_density(pft))*4/pi)**(pipe_tune3(pft)/2))*&
                  maxheight)**(1/(pipe_tune3(pft)/2+1))
             ! Vertical crown diameter derived from the height of that individual
             dia_ver(icir) = crown_to_height(pft)*(((woodmass_ind/(tree_ff(pft)*&
                  pipe_density(pft))*4/pi)**(pipe_tune3(pft)/2))*&
                  maxheight)**(1/(pipe_tune3(pft)/2+1))
             
          ENDIF

          ! Calculate the crown volume making use of the prescribed crown diameters
          ! Notice that we don't have to take the surface area of the stand
          ! into account.  total_bm will be in units of gC / m**2, while
          ! the canopy volume should be in units of m**3 / m**2. The surface
          ! area would cancel out when we take the ratio. The unit of total_cn_vol
          ! is thus similar to how precipitation is being expressed.
          ! Right now, total_cn_vol is in m**3/tree * tree/m**2 thus m
          cn_vol = pi/6*dia_hor(icir)*dia_hor(icir)*dia_ver(icir)
          total_cn_vol=total_cn_vol+cn_vol*ind(icir)  
       ENDDO

       ! Space available to be filled with the crown volumes. The same 
       ! unit convention as above is used (m)
       max_height = MAXVAL(wood_to_height(biomass_temp(:,:),pft,maxheight))
       

       IF (total_cn_vol.GT.max_height) THEN
       
          ! The prescribed crown diameters take too much space. Adjust the
          ! horizontal diameters such that the space is more correctly filled
          ! Given that we have to deal with different tree sizes and that the
          ! canopies are assumed to be ellipsoids this is a very complex problem
          ! to calculate exactly. We will make a couple of shortcuts.
          ! First we will assume that by packing ellipsoids we have the same
          ! efficiency as close-packing of equal spheres which is around 0.74.
          ! total_cn_vol = 1/6*pi*dia_hor^2*dia_ver
          ! dia_hor = crown_vertohor_dia*dia_ver
          ! <=> total_cn_vol = 1/6*pi*crown_vertohor_dia^2*dia_ver^3
          ! Reduce total_cn_vol to max_height * 0.74 (close packing) and calculate
          ! the adjusted crown_vertohor_dia.
          ! As an alternative scale both the vertyical and horizontal crown dimensions
          adjusted_crown_to_height = &
               ((max_height*crown_packing*6) / &
               (pi*(crown_vertohor_dia(pft)**2)*SUM(ind(:)*dia_ver(:)**3)))**(1./3.)
          dia_ver(:) = adjusted_crown_to_height * crown_to_height(pft) * &
               wood_to_height(biomass_temp(:,:),pft,maxheight)
          dia_hor(:) = crown_vertohor_dia(pft) * dia_ver(:)  

          ! Error checking
          total_cn_vol = zero
          DO icir = 1,ncirc
             total_cn_vol = total_cn_vol + &
                  pi/6*dia_hor(icir)*dia_hor(icir)*dia_ver(icir)*ind(icir)
          END DO
          IF (ABS(total_cn_vol-crown_packing*max_height).GT.min_stomate) THEN
             WRITE(numout,*) 'ind, ', ind(:)
             WRITE(numout,*) 'dia_ver, ', dia_ver(:)
             WRITE(numout,*) 'dia_hor, ', dia_hor(:)
             WRITE(numout,*) 'adjusted_crown_vertohor_dia, ',adjusted_crown_vertohor_dia
             WRITE(numout,*) 'adjusted_crown_to_height, ',adjusted_crown_to_height
             WRITE(numout,*) 'recalculated total_cn_vol, ',total_cn_vol
             WRITE(numout,*) 'crown_packing*max_height, ', crown_packing*max_height
             CALL ipslerr_p(3,'problem with crown packing',&
                  'adjusted horizonatal diameter is wrong',&
                  'check the calculations in the code','')
          END IF

       END IF

    ELSE

       WRITE(numout,*) 'ERROR: wood_to_height, pft ',pft
       CALL ipslerr_p (3,'wood_to_height', &
            'wood_to_height is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF  
   
 END SUBROUTINE wood_to_cn_dia


!! ================================================================================================================================
!! FUNCTION     : wood_to_cv
!!
!>\BRIEF        Calculate crown volume from woody biomass making use of the vertical and horizontal crown diameters
!!
!! DESCRIPTION : Calculate the volume of an elipsoid. V = 1/6 * pi * dia_hor**2 * dia_ver
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : crown volume (m3 ind-1)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  FUNCTION wood_to_cv(biomass_temp, ind, pft, maxheight)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std)                                          :: pft               !! PFT number (-)
    REAL(r_std), DIMENSION(:,:)                             :: biomass_temp      !! Biomass of an individual tree within a circ 
                                                                                 !! class @tex $(m^{2} ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:  )                             :: ind               !! Number of individuals within a circ class
                                                                                 !! @tex $(ind m^{2})$ @endtex
    REAL(r_std)                                             :: maxheight         !! pipe_tune2 for this pixel and pft
    
    !! 0.2 Output variables
                
    REAL(r_std), DIMENSION(ncirc)                           :: wood_to_cv        !! Crown volume of an individual tree 
                                                                                 !! @tex $(m^{2} ind{-1})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                          :: l                 !! Index
    REAL(r_std)                                             :: woodmass_ind      !! Woodmass of an individual tree
                                                                                 !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                           :: dia_ver           !! vertical crown diameter of an individual tree (m)
    REAL(r_std), DIMENSION(ncirc)                           :: dia_hor           !! horizontal crown diameter of an individual tree (m)

!_ ================================================================================================================================
 
 !! 1. Calculate crown volume from woodmass

    IF ( is_tree(pft) ) THEN

       ! Crown diameter of the individual tree (m2 ind-1)  
       CALL wood_to_cn_dia(biomass_temp, ind, pft, dia_hor, dia_ver,maxheight)
       wood_to_cv(:) = pi/6.*dia_ver(:)*dia_hor(:)*dia_hor(:)  
       ! WRITE(numout,*) 'wood_to_cv_eff, ', wood_to_cv_eff(:)

    ELSE

       WRITE(numout,*) 'pft ',pft
       CALL ipslerr_p (3,'wood_to_cv_eff', &
            'wood_to_cv_eff is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

  END FUNCTION wood_to_cv


!! ================================================================================================================================
!! FUNCTION     : Nmax
!!
!>\BRIEF        This function determines the maximum number of trees per hectare for 
!! a given quadratic mean diameter (Dg). It applies the self-thinning principle 
!! of Reineke (1933), with Dg instead of mean diameter (Dhote 1999).
!! Parameterization: Dhote (1999) for broad-leaved and Vacchiano (2008)
!! for needle-leaved.
!!
!! DESCRIPTION : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : Nmax
!!
!! REFERENCE(S)   : See above, module description.
!!
!! FLOWCHART   : None
!! \n
!_ ================================================================================================================================
   
FUNCTION Nmax(Dg,alpha_st,no_pft, ifm)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL(r_std)      :: Dg              !! Quadratic mean diameter (cm)
    INTEGER(i_std)   :: no_pft          !! Plant functional type (unitless) 
    REAL(r_std)      :: alpha_st        !! Alpha_self_thinning for specific ipts and ivm (unitless)
    INTEGER(i_std), INTENT(IN)         :: ifm  ! Forest management identifiers

    
    !! 0.2 Output variables

    REAL(r_std)      :: Nmax            !! Maximum number of trees according to the self-thinning model in (trees ha-1)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================
 !! 1. maximum number of trees per hectare for a given quadratic mean diameter

    IF (is_tree(no_pft)) THEN

       ! thinning curve of the MTC
       Nmax = (Dg/alpha_st)**(un/beta_self_thinning(no_pft))

       ! Truncate the range to avoid huge numbers due the exponental model used to describe self-thinning
       ! Don't set to ntreesmax else stomate_mark_kill will assume we are at self-thinning and
       ! won't apply background mortality. 
       ! Nmax = MIN( Nmax, ntreesmax(no_pft))
       Nmax = MAX( Nmax, dens_target(no_pft, ifm) )

    ELSE

       WRITE(numout,*) 'Self thinning is not defined for PFT, ', no_pft
       CALL ipslerr_p (3,'nmax', &
            'Self thinning is not defined for this PFT.', &
            'See the output file for more details.','')

    ENDIF

  END FUNCTION Nmax

!! ================================================================================================================================
!! FUNCTION     : calculate_rdi
!!
!>\BRIEF        Calculate the relative density index of a stand
!!
!! DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : calculate_rdi
!!
!! REFERENCE(S)   : Bellassen, V., Le Maire, G., Dhote, J.F., Viovy, N., Ciais, P., 2010. 
!! Modeling forest management within a global vegetation model – Part 1: 
!! model structure and general behaviour. Ecological Modelling 221, 2458–2474.
!!
!! FLOWCHART   : None
!! \n
!_ ================================================================================================================================
  
FUNCTION calculate_rdi(diameters_temp,circ_class_n_temp,alpha_st,ivm,ifm)

 !! 0. Variable and parameter declaration

    IMPLICIT NONE

    !! 0.1 Input variables
    REAL(r_std),DIMENSION(ncirc),INTENT(in)            :: circ_class_n_temp   !! The distribution of individuals
                                                                              !! for a given grid square and PFT  
    REAL(r_std),DIMENSION(ncirc),INTENT(in)            :: diameters_temp      !! The diameters of the trunks
                                                                              !! for each circ class (m)
    INTEGER(i_std),INTENT(IN)                          :: ivm,ifm             !! The number of the PFT
                                                                              !! (unitless)

    REAL(r_std),INTENT(IN)                             :: alpha_st            !! alpha_self_thinning for ipts x ivm (unitless)

    !! 0.2 Output variables
    REAL(r_std)                                        :: calculate_rdi       !! The RDI of the stand 
                                                                              !! (unitless, always positive)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                        :: Dg                  !! The quadratic mean diameter
    REAL(r_std)                                        :: ind_loc
    INTEGER(i_std)                                     :: icir                !! Indices

!_ ================================================================================================================================

 !! 1. Calculate the relative density index for a stand, based on
 !!    the number of trees in each circumference class and the
 !!    the diameter of a model tree in each class. 
    ind_loc=SUM(circ_class_n_temp(:))

    Dg=zero

    DO icir=1,ncirc
       Dg=Dg+circ_class_n_temp(icir)*m2_to_ha*(diameters_temp(icir)*m_to_cm)**2
    ENDDO
    Dg=SQRT(Dg/(ind_loc*m2_to_ha))
   
    calculate_rdi=(ind_loc*m2_to_ha)/Nmax(Dg,alpha_st,ivm,ifm)

  END FUNCTION calculate_rdi

!! ================================================================================================================================
!! FUNCTION     : weibull_class_dist
!!
!>\BRIEF        Define the main function to calculate the circonference class distribution
!               using weibull PDF and CDF
!!
!! DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : weibull_class_dist
!!
!! REFERENCE(S)   : -
!!
!! FLOWCHART   : None
!! \n
!_ ================================================================================================================================

!! ================================================================================================================================
!! FUNCTION     : cc_grid_max_x
!!
!>\BRIEF        Upper bound of the normalised circumference grid, per management regime
!!
!! DESCRIPTION  : The circumference classes sit on a grid whose multipliers are
!!                1/ncirc + (i-1)*max_x/ncirc, so the stand SPAN -- largest class diameter
!!                over smallest -- is 1+(ncirc-1)*max_x. With max_x hard-coded the span is
!!                tied to the number of classes, and refining ncirc silently describes a
!!
!! RETURN VALUE : cc_grid_max_x (unitless)
!! \n
!_ ================================================================================================================================

FUNCTION cc_grid_max_x(management)

      IMPLICIT NONE

      ! Input
      INTEGER(i_std), INTENT(in) :: management     !! Forest management type (unitless)
      ! Output
      REAL(r_std)                :: cc_grid_max_x  !! Upper bound of the normalised grid (unitless)
      ! Local
!_
!================================================================================================================================

   SELECT CASE (management)
      CASE (ifm_none, ifm_uneven)
         cc_grid_max_x = trois
      CASE (ifm_thin, ifm_cop, ifm_src)
         cc_grid_max_x = deux
      CASE DEFAULT
         CALL ipslerr_p(3,'cc_grid_max_x: unknown forest management system', &
              'the circumference grid bound is not defined for it', &
              'define your new management system in constantes_mtc', &
              'pft_parameters and pft_parameters_var')
   END SELECT


END FUNCTION cc_grid_max_x


FUNCTION weibull_class_dist(dia_lambda, ivm, management)

      IMPLICIT NONE

      ! Inputs  
      REAL, intent(in)                       :: dia_lambda          !! weibull shape parameter
      INTEGER(i_std),INTENT(IN)              :: ivm                 !! The number of the PFT
      INTEGER(i_std),INTENT(IN)              :: management          !! Forest management type (unitless)
      ! Outputs
      REAL, dimension(ncirc)                 :: weibull_class_dist  !! circumf class distribution using the weibull PDF

      ! Local
      INTEGER(i_std)                         :: i                   !! loop counter 
      REAL(r_std)                            :: dia,res             !! diameter increment, current diameter
      REAL(r_std), DIMENSION(ncirc)          :: CDF                 !! cumulative distribution function for circonf class 
      REAL(r_std)                            :: k, max_x            !! shape parameter and sum of CDF values
!_
!================================================================================================================================

   SELECT CASE (management)
      ! Calculate shape parameter k
      CASE (ifm_none) ! Unmanaged forest
         k = circ_dist_shape_none(ivm)
      CASE(ifm_thin) ! Only Even-age forest has an active diameter distribution management
         k = dia_lambda*circ_dist_shape_thin(ivm)/(largest_tree_dia(ivm)*100.0) + 0.5
      CASE(ifm_cop) ! short rotation plantation
         k = circ_dist_shape_cop(ivm)
      CASE(ifm_src) ! coppice forest
         k = circ_dist_shape_src(ivm)
      CASE(ifm_uneven) ! unevenage forest
         k = circ_dist_shape_uneven(ivm)
      CASE DEFAULT
        WRITE(numout,*) 'management is, ', management
         CALL ipslerr_p(3,'Shape parameter k is not known', &
              'for the forest management system', &
              'define your new management system in constantes_mtc', &
              'pft_parameters and pft_parameters_var')
   END SELECT

   ! Guillaume M. -- The grid bound now comes from cc_grid_max_x, single source shared with
   ! stomate_prescribe. The two used to carry their own copy of the same SELECT CASE, so a
   ! change to one silently described a different stand than the other.
   max_x = cc_grid_max_x(management)

   ! A higer value of max_x will increase the tail of the distribution.
   ! that's why we use a higer value of x-max for unmanaged and uneven forest 
   res = max_x / ncirc
   ! Initialize CDF and weibull_class_dist arrays to zero
   CDF = zero
   weibull_class_dist(:) = zero

   ! Calculate CDF values for each circonference class
   dia = 0.01
   DO i=1, size(CDF)
      dia = dia + res
      CDF(i) = (un - EXP(-(dia)**k))-(un - EXP(-(dia-res)**k))
   END DO

   DO i=1,ncirc
      ! rescale to ensure the sum of the CDF values equal to one.
      weibull_class_dist(i) = CDF(i)/SUM(CDF)
      IF (weibull_class_dist(i) <= zero) THEN
         WRITE(numout,*) 'weibull_class_dist for diameter class', i, &
              'is equal to zero or negative'
         CALL ipslerr_p(3,'weibull_class_dist is equal to zero or negative',&
              'this should not be the case', 'check what went wrong', '')
      END IF
   END DO
END FUNCTION



!
!================================================================================================================================
!! FUNCTION     : calculate_rdi_boundaries
!!
!>\BRIEF        calculate upper and lower boundaries for some resource demand
!               index (RDI) based on input parameters
!!
!! DESCRIPTION   : The subroutine takes several input arguments:
!!                  dia: A real number representing diameter.
!!                  ivm and ifm: Integer values representing PFT (Plant Functional Type) and forest
!!                               management.
!!                  upper and lower: Real numbers passed by reference (inout) to store the upper and
!!                                   lower boundaries. It initializes the upper and lower variables to zero.
!!
!!                It calculates the value of d by subtracting 1 from the size of poly_rdi(ivm,
!!                  ifm, 1,:), where poly_rdi is presumably a multi-dimensional array. This suggests
!!                  that poly_rdi contains polynomial coefficients for some calculation related to
!!                  the RDI.
!!
!!                It enters a loop that runs from 0 to d. Inside the loop, it calculates the upper
!!                  and lower boundaries using polynomial coefficients from poly_rdi and the value
!!                  of diameters
!!
!!                After the loop, it sets the upper and lower values to the maximum of themselves
!!                  and zero. This ensures that both upper and lower values are non-negative.
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : upper and lower rdi boundaries 
!!
!! REFERENCE(S)   : -
!!
!! FLOWCHART   : None
!! \n
!_
!================================================================================================================================
SUBROUTINE calculate_rdi_boundaries(dia, ivm, ifm, upper, lower)
  REAL, intent(in)                   :: dia     ! Input diameter in meters
  INTEGER(i_std), INTENT(IN)         :: ivm, ifm  ! PFT and forest management identifiers
  REAL, intent(inout)                :: upper, lower  ! Output upper and lower RDI boundaries
  INTEGER(i_std)                     :: i, d     ! Loop counters and dimensionality variable

  upper = 0.0  ! Initialize upper boundary
  lower = 0.0  ! Initialize lower boundary

  d = SIZE(poly_rdi(ivm, ifm, 1, :)) - 1  ! Calculate the degree of the polynomial

  ! Check if the diameter is in meters and not centimeters
  if (dia > deux) then
    WRITE(numout,*) 'For PFT ', ivm, ' at forest management ', ifm
    CALL ipslerr_p(3, 'dia must be in meters', &
      'It appears you have already converted it to cm.', &
      'Please change the unit of dia from cm to m.', '')
  end if

  ! Calculate upper and lower RDI boundaries using polynomial coefficients
  DO i = 0, d
    upper = upper + poly_rdi(ivm, ifm, 1, i+1) * (dia * m_to_cm) ** i
    lower = lower + poly_rdi(ivm, ifm, 2, i+1) * (dia * m_to_cm) ** i
  END DO

  ! Ensure that boundaries are non-negative
  upper = MAX(upper, 0.0)
  lower = MAX(lower, 0.0)

  ! Check if upper boundary is lower than lower boundary (an error condition)
  if (upper < lower) then
    WRITE(numout,*) 'For PFT ', ivm, ' at forest management ', ifm
    WRITE(numout,*) 'Upper ', upper
    WRITE(numout,*) 'Lower ', lower
    WRITE(numout,*) 'dia ', dia
    WRITE(numout,*) 'dimension ', d
    CALL ipslerr_p(3, 'rdi_target_upper is lower than rdi_target_lower', &
      'This should not be the case. Check what went wrong.', '','')
  end if

END SUBROUTINE



!! ================================================================================================================================
!! SUBROUTINE  : check_vegetation_area
!!
!>\BRIEF       
!!
!! DESCRIPTION : 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_vegetation_area(check_point, npts, veget_max_begin, &
       veget_max, test_type, change_nobio)

  !! 0. Variable and parameter description

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts                 !! Domain size (unitless)
    CHARACTER(*),INTENT(in)                           :: check_point          !! A flag to indicate at which
                                                                              !! point in the code we're doing
                                                                              !! this check
    CHARACTER(*), INTENT(in)                          :: test_type            !! Which test do we want to run on veget_max
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max            !! "maximal" coverage fraction of a PFT 
                                                                              !! (LAI -> infinity) on ground. May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless; 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max_begin      !! temporary storage of veget_max to check area conservation
    REAL(r_std), DIMENSION(:), INTENT(in), OPTIONAL   :: change_nobio         !! Change in the nobio fraction due to urbanisation or possibly glacier retreat                

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: iele, ipts, ivm      !! Indices
    INTEGER(i_std)                                    :: ipar, icir           !! Indices
    INTEGER(i_std)                                    :: ivma,iagec,ipft     !! Counter
    REAL(r_std)                                       :: temp_begin, temp_end !! temporary variables to accumulate veget_max of
                                                                              !! different age classes of the same species  

!_ ================================================================================================================================

    SELECT CASE (test_type)
       
    CASE ('pft')
       
       ! In all subroutines where veget_max should be preserved at the
       ! PFT level, the most strict test can be performed
       DO ipts=1,npts

          DO ivm = 1,nvm

             IF (ABS(veget_max_begin(ipts,ivm) - veget_max(ipts,ivm)) .GT. &
                  min_stomate ) THEN
                WRITE(numout,*) 'Surface area is not preserved in ', &
                     TRIM(check_point)
                WRITE(numout,*) 'Area has changed for PFT, ',ipts, ivm
                WRITE(numout,*) 'Change from, ',veget_max_begin(ipts,ivm), &
                     ' to ', veget_max(ipts,ivm)
                IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                        'Surface area not preserved','','')
             ENDIF
             
          ENDDO !nvm

          ! One message per pixel
!!$          IF (err_act.GT.1) WRITE(numout,*) 'Surface area preserved in ',  &
!!$               TRIM(check_point), ipts
             
       ENDDO ! npts

    CASE ('ageclass')

       ! In subroutines where veget_max can move from one pft to 
       ! another within the same species. This is the test to use
       DO ipts=1,npts
          
          DO ivma=1,nvmap

             ! get pft number
             ivm=start_index(ivma)

             ! If we only have a single age class for this
             ! PFT, we can simplify the test
             IF(nagec_pft(ivma) .EQ. 1) THEN

                IF (veget_max_begin(ipts,ivm) - veget_max(ipts,ivm) .GT. &
                     min_stomate ) THEN
                   WRITE(numout,*) 'Surface area is not preserved in ', &
                        TRIM(check_point)
                   WRITE(numout,*) 'Area has changed for PFT, ',ipts, ivm
                   WRITE(numout,*) 'Change from, ',veget_max_begin(ipts,ivm), &
                        ' to ', veget_max(ipts,ivm)
                   IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                        'Surface area not preserved','','')
                ENDIF

             ELSE ! more than 1 age class

                temp_begin = 0
                temp_end = 0
                DO iagec = nagec_pft(ivma),1,-1
                   
                   ! Accumulate veget_max for all pfts/age classes
                   ! that belong to the same species
                   ipft = ivm+iagec-1
                   temp_begin = temp_begin + veget_max_begin(ipts,ipft)
                   temp_end = temp_end + veget_max(ipts,ipft)

                ENDDO
                
                IF (ABS(temp_begin - temp_end) .GT. min_stomate ) THEN
                   WRITE(numout,*) 'Surface area is not preserved in ', &
                        TRIM(check_point)
                   WRITE(numout,*) 'Area has changed for PFT, ',ipts, ivm
                   WRITE(numout,*) 'Change from, ',veget_max_begin(ipts,ivm), &
                        ' to ', veget_max(ipts,ivm)
                   IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                        'Surface area not preserved','','')
                ENDIF

             ENDIF ! number of age classes

          ENDDO ! nvmap

          ! One message per pixel
          !IF (err_act.GT.1) WRITE(numout,*) 'Surface area preserved in ',  &
          !     TRIM(check_point), ipts
          
       ENDDO ! ipts

    CASE('pixel')

       ! Use this test in subroutines where veget_max can be moved
       ! between species. This is the least strict test and only 
       ! test wheter the total veget_max was preserved.
       DO ipts=1,npts
          
          ! Change_nobio is the change in the frac_nobio that may occur in
          ! LCC. This change should be accounted for in the mass balance check.
          ! ::change_nobio is defined as an optional element and is only present
          ! when the surface area is checked in sapiens_lcchange
          IF(PRESENT(change_nobio))THEN
             IF (ABS(SUM(veget_max_begin(ipts,:)) - SUM(veget_max(ipts,:)) - change_nobio(ipts)) .GT. &
                  min_stomate ) THEN
                WRITE(numout,*) 'Surface area is not preserved in ', &
                     TRIM(check_point)
                WRITE(numout,*) 'Area in pixel ', ipts,' has changed form, ', & 
                     SUM(veget_max_begin(ipts,:)), 'to,', SUM(veget_max(ipts,:))
                IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                     'Surface area not preserved','','')
             ENDIF
          ELSE
             IF (ABS(SUM(veget_max_begin(ipts,:)) - SUM(veget_max(ipts,:))) .GT. &
                  min_stomate ) THEN
                WRITE(numout,*) 'Surface area is not preserved in ', &
                     TRIM(check_point)
                WRITE(numout,*) 'Area in pixel ', ipts,' has changed form, ', & 
                     SUM(veget_max_begin(ipts,:)), 'to,', SUM(veget_max(ipts,:))
                IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                     'Surface area not preserved','','')
             ENDIF
          ENDIF

       ENDDO ! npts
       
    CASE DEFAULT

        WRITE(numout,*) 'Error: an unknown type/resolution was requested'
        WRITE(numout,*) ' for the surface area check. Check the cases'
        WRITE(numout,*) ' defined in src_parameters/function_library or'
        WRITE(numout,*) ' the spelling of the type in the CALL argument',&
            test_type
        CALL ipslerr_p (3,'surface_area_check in function_library',&
             'unknown type', 'Check *_out_orchidee for more details','')        
       
    END SELECT

  END SUBROUTINE check_vegetation_area


!! ================================================================================================================================
!! SUBROUTINE  : check_area_invariant
!!
!>\BRIEF        Assert, SLOT BY SLOT, that a slot either carries area or is exactly empty.
!!
!! DESCRIPTION  : The existing area guards compare SUMS over a pixel. An area moved between two
!! slots of the same group is null in that projection, so it is invisible to them. This routine
!! looks at every slot instead, and reports two violations:
!!   (a) 0 < veget_max < min_stomate      : a residue left by a transfer meant to empty the slot
!!   (b) veget_max <= min_stomate and the slot still carries individuals: the state that stops
!!       stomate_prescribe on "Unexpected condition" the following year
!! It also counts the slots that fall inside the window check_veget claims to reject, because it
!! is not established that check_veget ever sees them.
!!
!! Diagnostic only: reads, WRITE and ipslerr_p(2). STRICTLY BIT-NEUTRAL.
!!
!! RECENT CHANGE(S) : created 2026-08-06, step 0 of design/MODULE_DESIGN_AREA_TRANSFER_AUDIT.md
!!
!! MAIN OUTPUT VARIABLE(S) : none
!!
!! REFERENCES   : None
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE check_area_invariant(check_point, npts, veget_max, circ_class_n)

  !! 0. Variable and parameter description

    !! 0.1 Input variables
    CHARACTER(*), INTENT(in)                          :: check_point          !! Where in the code this check runs
    INTEGER(i_std), INTENT(in)                        :: npts                 !! Domain size (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max            !! Maximal coverage fraction of a PFT
                                                                              !! (unitless; 0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: circ_class_n         !! Number of individuals per circumference
                                                                              !! class @tex $(ind m^{-2})$ @endtex

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: ipts, ivm            !! Indices
    INTEGER(i_std)                                    :: n_residue            !! Slots with area below min_stomate
    INTEGER(i_std)                                    :: n_orphan             !! Slots with content but no area
    INTEGER(i_std)                                    :: n_window             !! Slots inside the check_veget window
    INTEGER(i_std)                                    :: n_shown              !! Detailed lines already written
    INTEGER(i_std), PARAMETER                         :: n_show_max = 10      !! Cap on detailed lines per call
    REAL(r_std)                                       :: dens                 !! Individuals carried by the slot

!_ ================================================================================================================================

    n_residue = 0
    n_orphan  = 0
    n_window  = 0
    n_shown   = 0

    DO ivm = 1,nvm
       DO ipts = 1,npts

          dens = SUM(circ_class_n(ipts,ivm,:))

          ! (a) An area below min_stomate is not an area. It is what a transfer meant to empty
          ! the slot left behind.
          IF (veget_max(ipts,ivm) .GT. zero .AND. veget_max(ipts,ivm) .LT. min_stomate) THEN
             n_residue = n_residue + 1
             IF (n_shown .LT. n_show_max) THEN
                WRITE(numout,*) '[AREA_INV] ', TRIM(check_point), ' residue  ipts', ipts, &
                     ' ivm', ivm, ' veget_max', veget_max(ipts,ivm), ' n', dens
                n_shown = n_shown + 1
             ENDIF
          ENDIF

          ! (b) Content without ground. THIS is the state that stops stomate_prescribe.
          IF (veget_max(ipts,ivm) .LE. min_stomate .AND. dens .GT. min_stomate) THEN
             n_orphan = n_orphan + 1
             IF (n_shown .LT. n_show_max) THEN
                WRITE(numout,*) '[AREA_INV] ', TRIM(check_point), ' orphan   ipts', ipts, &
                     ' ivm', ivm, ' veget_max', veget_max(ipts,ivm), ' n', dens
                n_shown = n_shown + 1
             ENDIF
          ENDIF

          ! Q2 of the audit: check_veget rejects this window, yet the measured residues sat in it
          ! without stopping the model. Count them here so we learn whether check_veget ever sees
          ! them at all.
          IF (veget_max(ipts,ivm) .GT. EPSILON(un)*100._r_std .AND. &
               veget_max(ipts,ivm) .LT. min_vegfrac) n_window = n_window + 1

       ENDDO
    ENDDO

    IF (n_residue .GT. 0 .OR. n_orphan .GT. 0) THEN
       WRITE(numout,*) '[AREA_INV] ', TRIM(check_point), ' TOTAL residue', n_residue, &
            ' orphan', n_orphan, ' in-check_veget-window', n_window
       CALL flush(numout)
       CALL ipslerr_p(2,'check_area_invariant', &
            'a slot carries area below min_stomate, or content with no area', &
            TRIM(check_point),'')
    ELSEIF (n_window .GT. 0) THEN
       WRITE(numout,*) '[AREA_INV] ', TRIM(check_point), ' in-check_veget-window', n_window
       CALL flush(numout)
    ENDIF

  END SUBROUTINE check_area_invariant


!! ================================================================================================================================
!! SUBROUTINE  : check_vegetation_area_point
!!
!>\BRIEF: check changes in veget_max for a single point: ipts x ivm      
!!
!! DESCRIPTION : 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_vegetation_area_point(check_point, veget_max_begin, veget_max)

  !! 0. Variable and parameter description

    !! 0.1 Input variables
    CHARACTER(*), INTENT(in)                          :: check_point          !! A flag to indicate at which
                                                                              !! point in the code we're doing
                                                                              !! this check
    REAL(r_std), INTENT(in)                           :: veget_max            !! "maximal" coverage fraction of a PFT 
                                                                              !! (LAI -> infinity) on ground. May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless; 0-1)
    REAL(r_std), INTENT(in)                           :: veget_max_begin      !! temporary storage of veget_max to check area conservation


    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: iele, ipts, ivm      !! Indices
    INTEGER(i_std)                                    :: ipar, icir           !! Indices
    INTEGER(i_std)                                    :: ivma,iagec,ipft     !! Counter
    REAL(r_std)                                       :: temp_begin, temp_end !! temporary variables to accumulate veget_max of
                                                                              !! different age classes of the same species  

!_ ================================================================================================================================


    ! Use this test in subroutines where the calculations are made over a single
    ! point: ipts x ivm
    IF (ABS(veget_max_begin - veget_max) .GT. min_stomate ) THEN
       WRITE(numout,*) 'Surface area is not preserved in ', &
            TRIM(check_point)
       WRITE(numout,*) 'Area has changed form, ', & 
            veget_max_begin, 'to,', veget_max
       IF(err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
            'Surface area not preserved','','')
    ENDIF
  
  END SUBROUTINE check_vegetation_area_point

!! ================================================================================================================================
!! SUBROUTINE  : check_mass balance
!!
!>\BRIEF       
!!
!! DESCRIPTION : 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_mass_balance(check_point, closure_intern, npts, pool_end, &
            pool_start, veget_max, test_type, pft)

 !! 0. Variable and parameter description

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: npts                 !! Domain size (unitless)
    CHARACTER(*),INTENT(in)                           :: check_point          !! A flag to indicate at which
                                                                              !! point in the code we're doing
                                                                              !! this check
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: closure_intern       !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: pool_start           !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: pool_end             !! Start and end pool of this routine 
    CHARACTER(*), INTENT(in)                          :: test_type            !! Which test do we want to run on veget_max
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max            !! "maximal" coverage fraction of a PFT 
                                                                              !! (LAI -> infinity) on ground. May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless; 0-1) 
    INTEGER(i_std), OPTIONAL, INTENT(in)               :: pft                 !! Number of the PFT for which the mass balance has to be checked 

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                           :: iele, ipts, ivm      !! Indices
    INTEGER                                           :: ipar, icir, err      !! Indices
    INTEGER(i_std)                                    :: ivma, iagec, ipft    !! Counter
    REAL(r_std)                                       :: temp_begin, temp_end !! temporary variables to accumulate veget_max of

    CHARACTER(LEN=10)                                 :: var_name             !! Name of element
    LOGICAL                                           :: ltemp                !! temporary logical variable
!_ ================================================================================================================================

    SELECT CASE (test_type)
       
    CASE ('ipts')

       ! Close the mass balance for one specific pixel and pft. 
       ! This :test is used for subroutines that are called within 
       ! a DO ipts = 1,npts and DO ivm = 1,nvm loop. The mass
       ! balance is checked for a single pft at a single pixel. 

       DO iele=1,nelements
          
          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'
          err = zero

          ! Mass balance closure      
          DO ipts=1,npts
             IF(ABS(SUM(closure_intern(ipts,:,iele))) .GT. min_stomate)THEN

                ! Count the number of errors
                err = err + 1

                ! Give extra details of where the mass balance
                ! problem occurs. Stop the model if the user
                ! asked for it
                IF(err_act.GT.1) THEN
                   WRITE(numout,*) 'Error: mass balance is not closed in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) ' Problem occurs in pixel, ',ipts
                   WRITE(numout,*) ' Difference is, ', &
                        SUM(closure_intern(ipts,:,iele))
                   WRITE(numout,*) ' pool_end,pool_start: ', &
                        SUM(pool_end(ipts,:,iele)), &
                        SUM(pool_start(ipts,:,iele))
                   CALL flush(numout)
     
                   ! Exit or write warning depening on plev
                   CALL ipslerr_p(plev, TRIM(check_point), &
                        'Mass balance error','','')
                ENDIF

             ENDIF
          END DO

          !+++CHECK+++
          ! A lot of output - only write a message if their is a problem
!!$          ! One message per pixel. The err condition is needed because under
!!$          ! err_act = 2 the model will write an error to the out_orchidee file
!!$          ! but will not be stopped. Therefore it should be checked whether we
!!$          ! have an error or not
!!$          IF (err_act.GT.1 .AND. err.EQ.0) THEN
!!$             WRITE(numout,*) 'Mass balance closure in ', &
!!$                  TRIM(check_point), ' for ',var_name
!!$          ENDIF
          !++++++++++

       ENDDO ! nelements

    CASE ('pft')
    
       ! Close mass balance at the pixel and pft level.
       DO iele=1,nelements

          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'

          DO ipts=1,npts

             err = zero

             DO ivm=1,nvm

                !check for NaN
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_start(ipts,ivm,iele))
#else
                ltemp=isnan(pool_start(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*)'At pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' pool_start is NaN'
                   WRITE(numout,*) ' pool_start: ', pool_start(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_end(ipts,ivm,iele))
#else
                ltemp=isnan(pool_end(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' pool_end is NaN'
                   WRITE(numout,*) ' pool_end: ', pool_end(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(closure_intern(ipts,ivm,iele))
#else
                ltemp=isnan(closure_intern(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' closure_intern is NaN'
                   WRITE(numout,*) ' closure_intern: ', closure_intern(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point), 'NaN problem','','')
                ENDIF

                ! Mass balance closure          
                IF(ABS(closure_intern(ipts,ivm,iele)) .GT. min_stomate)THEN
              
                   err = err+1
      
                   ! If the model is under development or being tested we will
                   ! output extra information to facilitate debugging.
                   IF(err_act.GT.1)THEN
                      WRITE(numout,*) 'Error: mass balance is not closed in ', &
                           TRIM(check_point), ' for ', var_name
                      WRITE(numout,*) ' ipts, ivm; ', ipts,ivm
                      WRITE(numout,*) ' Difference is, ', &
                           closure_intern(ipts,ivm,iele)
                      WRITE(numout,*) ' pool_end,pool_start: ', &
                           pool_end(ipts,ivm,iele), &
                           pool_start(ipts,ivm,iele)
                      CALL flush(numout)
                      ! Exit or write warning depening on plev
                      CALL ipslerr_p(plev, TRIM(check_point), &
                           'Mass balance error','','')
                   ENDIF
                ENDIF
             ENDDO ! nvm

             !+++CHECK+++
             ! A lot of output - only write a message if their is a problem
!!$             ! One message per pixel
!!$             IF (err_act.GT.1 .AND. err.EQ.zero) THEN
!!$                WRITE(numout,*) 'Mass balance closure in ', &
!!$                     TRIM(check_point), ' for ',var_name
!!$             ENDIF
             !++++++++++

          ENDDO ! npts

       ENDDO ! nelements

    CASE ('ipft')
    
       ! Close mass balance at the pixel and pft level. Check for a
       ! specific pft. This option is called from within a DO-loop.
       DO iele=1,nelements

          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'

          DO ipts=1,npts
             
             err = zero

             !check for NaN
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_start(ipts,pft,iele))
#else
                ltemp=isnan(pool_start(ipts,pft,iele))
#endif
             IF (ltemp) THEN
                WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                WRITE(numout,*)'At the pixel',ipts,' for element',iele,&
                     ' and for PFT',pft,' pool_start is NaN'
                WRITE(numout,*) ' pool_start: ', pool_start(ipts,:,iele)
                CALL flush(numout)
                CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
             ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_end(ipts,pft,iele))
#else
                ltemp=isnan(pool_end(ipts,pft,iele))
#endif
             IF (ltemp) THEN
                WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                     ' and for PFT',pft,' pool_end is NaN'
                WRITE(numout,*) ' pool_end: ', pool_end(ipts,:,iele)
                CALL flush(numout)
                CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
             ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(closure_intern(ipts,pft,iele))
#else
                ltemp=isnan(closure_intern(ipts,pft,iele))
#endif
             IF (ltemp) THEN
                WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                     ' and for PFT',pft,' closure_intern is NaN'
                WRITE(numout,*) ' closure_intern: ', closure_intern(ipts,:,iele)
                CALL flush(numout)
                CALL ipslerr_p(plev,TRIM(check_point), 'NaN problem','','')
             ENDIF
             
             
             ! Mass balance closure          
             IF(ABS(closure_intern(ipts,pft,iele)) .GT. min_stomate)THEN
                
                err = err+1
                
                ! If the model is under development or being tested we will
                ! output extra information to facilitate debugging.
                IF(err_act.GT.1)THEN
                   WRITE(numout,*) 'Error: mass balance is not closed in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) ' ipts, pft; ', ipts,pft
                   WRITE(numout,*) ' Difference is, ', &
                        closure_intern(ipts,pft,iele)
                   WRITE(numout,*) ' pool_end,pool_start: ', &
                        pool_end(ipts,pft,iele), &
                        pool_start(ipts,pft,iele)
                   CALL flush(numout)
                   ! Exit or write warning depening on plev
                   CALL ipslerr_p(plev, TRIM(check_point), &
                        'Mass balance error','','')
                ENDIF
             ENDIF

          ENDDO ! npts

       ENDDO ! nelements

    CASE ('ageclass')
       
       ! This option is only used when err_act.GE.2

       ! In the subroutine under test, biomass can move
       ! from one age class to another age class in the
       ! same species. Test across all age classes for
       ! the same species.
       DO iele=1,nelements

          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'
        
          DO ipts=1,npts

             err = zero

             DO ivma=1,nvmap

                ! get pft number
                ivm=start_index(ivma)

                !check for NaN
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_start(ipts,ivm,iele))
#else
                ltemp=isnan(pool_start(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*)'At the pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' pool_start is NaN'
                   WRITE(numout,*) ' pool_start: ', pool_start(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(pool_end(ipts,ivm,iele))
#else
                ltemp=isnan(pool_end(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' pool_end is NaN'
                   WRITE(numout,*) ' pool_end: ', pool_end(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(closure_intern(ipts,ivm,iele))
#else
                ltemp=isnan(closure_intern(ipts,ivm,iele))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' and for PFT',ivm,' closure_intern is NaN'
                   WRITE(numout,*) ' closure_intern: ', closure_intern(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point), 'NaN problem','','')
                ENDIF

                ! If we only have a single age class for this
                ! PFT, we can simplify the test
                IF(nagec_pft(ivma) .EQ. 1) THEN

                   ! Mass balance closure          
                   IF(ABS(closure_intern(ipts,ivm,iele)) .GT. min_stomate)THEN
                      err=err+1
                      WRITE(numout,*) 'Error: mass balance is not closed in ', &
                           TRIM(check_point), ' for ', var_name
                      WRITE(numout,*) ' ipts, ivm; ', ipts,ivm
                      WRITE(numout,*) ' Difference is, ', &
                           closure_intern(ipts,ivm,iele)
                      WRITE(numout,*) ' pool_end,pool_start: ', &
                           pool_end(ipts,ivm,iele), &
                           pool_start(ipts,ivm,iele)
                      CALL flush(numout)
                      IF (err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                           'Mass balance error','','')

                   ENDIF ! test for closure

                ELSE ! more than 1 age class

                   temp_begin = 0
                   DO iagec = nagec_pft(ivma),1,-1
                      
                      ! Accumulate check_intern for all pfts/age classes
                      ! that belong to the same species
                      ipft = ivm+iagec-1
                      temp_begin = temp_begin + closure_intern(ipts,ipft,iele)
                      
                   ENDDO ! nagec_pft
                   
                   IF (ABS(temp_begin) .GT. min_stomate) THEN

                      ! mass balance
                      err=err+1
                      WRITE(numout,*) 'Error: mass balance is not closed in ', &
                           TRIM(check_point), ' for ', var_name
                      WRITE(numout,*) ' ipts, ivma; ', ipts,ivma
                      WRITE(numout,*) ' first and last pft, ', &
                           start_index(ivma), ipft
                      WRITE(numout,*) ' Difference is, ', &
                           temp_begin
                      CALL flush(numout)
                      IF (err_act.GT.1) CALL ipslerr_p (plev, TRIM(check_point), &
                           'Mass balance error','','')

                   ENDIF ! test for closure

                ENDIF ! number of age classes

             ENDDO ! ivmap

             
             !+++CHECK+++
             ! A lot of output - only write a message if their is a problem
!!$             ! One message per pixel
!!$             IF (err_act.GT.1 .AND. err.EQ.zero) THEN
!!$                WRITE(numout,*) 'Mass balance closure in ', &
!!$                     TRIM(check_point), ' for ',var_name
!!$             ENDIF
             !++++++++++

          ENDDO ! npts

       ENDDO ! iele

    CASE ('pixel')

       ! For testing subroutines where biomass can also move
       ! from one species to another.
       DO iele = 1, nelements

          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'
        
          DO ipts=1,npts

             err = zero
         
                !check for NaN
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(SUM(pool_start(ipts,:,iele)))
#else
                ltemp=isnan(SUM(pool_start(ipts,:,iele)))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*)'At the pixel',ipts,' for element',iele,&
                        ' pool_start is NaN'
                   WRITE(numout,*) ' pool_start: ', pool_start(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(SUM(pool_end(ipts,:,iele)))
#else
                ltemp=isnan(SUM(pool_end(ipts,:,iele)))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' pool_end is NaN'
                   WRITE(numout,*) ' pool_end: ', pool_end(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
                ENDIF
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(SUM(closure_intern(ipts,:,iele)))
#else
                ltemp=isnan(SUM(closure_intern(ipts,:,iele)))
#endif
                IF (ltemp) THEN
                   WRITE(numout,*) 'Error: found an NaN in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) 'At the pixel',ipts,' for element',iele,&
                        ' closure_intern is NaN'
                   WRITE(numout,*) ' closure_intern: ', closure_intern(ipts,:,iele)
                   CALL flush(numout)
                   CALL ipslerr_p(plev,TRIM(check_point), 'NaN problem','','')
                ENDIF
                
                ! Mass balance closure          
                IF(ABS(SUM(closure_intern(ipts,:,iele))) .GT. min_stomate)THEN

                err = err+1

                IF(err_act.GT.1)THEN
                   WRITE(numout,*) 'Error: mass balance is not closed in ', &
                        TRIM(check_point), ' for ', var_name
                   WRITE(numout,*) ' ipts, ', ipts
                   WRITE(numout,*) ' Difference is, ', &
                        SUM(closure_intern(ipts,:,iele))
                   WRITE(numout,*) ' pool_end,pool_start: ', &
                        SUM(pool_end(ipts,:,iele)), &
                        SUM(pool_start(ipts,:,iele))
                   CALL flush(numout)
                   ! Exit or write warning depening on plev
                   CALL ipslerr_p(plev, TRIM(check_point), &
                        'Mass balance error','','')
                ENDIF
                
             ENDIF

             !+++CHECK+++
             ! A lot of output - only write a message if their is a problem
!!$             ! One message per pixel
!!$             IF (err_act.GT.1 .AND. err.EQ.zero) THEN
!!$                WRITE(numout,*) 'Mass balance closure in ', &
!!$                     TRIM(check_point), ' for ',var_name
!!$             ENDIF
             !++++++++++
             
          ENDDO ! npts
          
       ENDDO ! iele

    CASE ('products')

       ! This option is only used when err_act.GE.2
      
       ! testing mass balance for the product pool. Veget_max is not used
       ! ivm is not defined for closure_intern 
       DO iele=1,nelements

          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'
          
          DO ipts=1,npts

             err = zero

             IF(ABS(closure_intern(ipts,1,iele)) .GT. min_stomate)THEN

                err = err+1 
                WRITE(numout,*) 'Error: mass balance is not closed'//&
                     ' in product use ', var_name
                WRITE(numout,*) 'ipts,ivm,iele ', ipts
                WRITE(numout,*) 'Difference is, ', &
                     closure_intern(ipts,1,iele)
                WRITE(numout,*) 'pool_end, pool_start: ', &
                     pool_end(ipts,1,iele), pool_start(ipts,1,iele)
                CALL flush(numout)
                IF(err_act.GT.1) CALL ipslerr_p (plev,'product use', &
                     'Mass balance error','','')
             ENDIF
             
             !+++CHECK+++
             ! A lot of output - only write a message if their is a problem
!!$             ! One message per pixel
!!$             IF (err_act.GT.1 .AND. err.EQ.zero) THEN
!!$                WRITE(numout,*) 'Mass balance closure in ', &
!!$                     TRIM(check_point), ' for ',var_name
!!$             ENDIF
             !++++++++++  

          ENDDO
       ENDDO
    
    CASE DEFAULT

        WRITE(numout,*) 'Error: an unknown type/resolution was requested'
        WRITE(numout,*) ' for the mass balance check. Check the cases'
        WRITE(numout,*) ' defined in src_parameters/function_library or'
        WRITE(numout,*) ' the spelling of the type in the CALL argument',&
             test_type
        CALL flush(numout)
        CALL ipslerr_p (3,'mass_balance_check in function_library',&
        'unknown type', 'Check *_out_orchidee for more details','')        
       
    END SELECT

 END SUBROUTINE check_mass_balance

 

!! ================================================================================================================================
!! SUBROUTINE  : check_mass balance_point
!!
!>\BRIEF       
!!
!! DESCRIPTION : 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE check_mass_balance_point(check_point, closure_intern, pool_end, &
            pool_start, veget_max)

 !! 0. Variable and parameter description

    !! 0.1 Input variables
    
    CHARACTER(*),INTENT(in)                           :: check_point          !! A flag to indicate at which
                                                                              !! point in the code we're doing
                                                                              !! this check
    REAL(r_std), DIMENSION(:), INTENT(in)             :: closure_intern       !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)             :: pool_start           !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)             :: pool_end             !! Start and end pool of this routine 
     REAL(r_std), INTENT(in)                          :: veget_max            !! "maximal" coverage fraction of a PFT 
                                                                              !! (LAI -> infinity) on ground. May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless; 0-1) 
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                           :: iele, ipts, ivm      !! Indices
    INTEGER                                           :: ipar, icir, err      !! Indices
    INTEGER(i_std)                                    :: ivma, iagec, ipft    !! Counter
    REAL(r_std)                                       :: temp_begin, temp_end !! temporary variables to accumulate veget_max of

    CHARACTER(LEN=10)                                 :: var_name             !! Name of element
    LOGICAL                                           :: ltemp                !! temporary logical variable
!_ ================================================================================================================================

       ! Close the mass balance for one specific pixel and pft. 
       ! This :test is used for subroutines that are called within 
       ! a DO ipts = 1,npts and DO ivm = 1,nvm loop. The mass
       ! balance is checked for a single pft at a single pixel. 
       DO iele=1,nelements
          
          IF (iele.EQ.icarbon) var_name = 'carbon'
          IF (iele.EQ.initrogen) var_name = 'nitrogen'
          err = zero

          !check for NaN
#ifdef CPP_IEEE_ARITHMETIC
          ltemp=IEEE_IS_NAN(pool_start(iele))
#else
          ltemp=isnan(pool_start(iele))
#endif
          IF (ltemp) THEN
             WRITE(numout,*) 'Checking Nan of a single ipts x ivm'
             WRITE(numout,*) 'Error: found an NaN in ', &
                  TRIM(check_point), ' for ', var_name
             WRITE(numout,*)' For element',iele,&
                  ' pool_start is NaN'
             WRITE(numout,*) ' pool_start: ', pool_start(iele)
             CALL flush(numout)
             CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
          ENDIF
#ifdef CPP_IEEE_ARITHMETIC
          ltemp=IEEE_IS_NAN(pool_end(iele))
#else
          ltemp=isnan(pool_end(iele))
#endif
          IF (ltemp) THEN
             WRITE(numout,*) 'Checking Nan of a single ipts x ivm'
             WRITE(numout,*) 'Error: found an NaN in ', &
                  TRIM(check_point), ' for ', var_name
             WRITE(numout,*) 'For element',iele,&
                  ' pool_end is NaN'
             WRITE(numout,*) ' pool_end: ', pool_end(iele)
             CALL flush(numout)
             CALL ipslerr_p(plev,TRIM(check_point),'NaN problem','','')
          ENDIF
#ifdef CPP_IEEE_ARITHMETIC
          ltemp=IEEE_IS_NAN(closure_intern(iele))
#else
          ltemp=isnan(closure_intern(iele))
#endif
          IF (ltemp) THEN
             WRITE(numout,*) 'Checking Nan of a single ipts x ivm'
             WRITE(numout,*) 'Error: found an NaN in ', &
                  TRIM(check_point), ' for ', var_name
             WRITE(numout,*) 'Ror element',iele,&
                  ' closure_intern is NaN'
             WRITE(numout,*) ' closure_intern: ', closure_intern(iele)
             CALL flush(numout)
             CALL ipslerr_p(plev,TRIM(check_point), 'NaN problem','','')
          ENDIF

          ! Mass balance closure      
          IF(ABS(closure_intern(iele)) .GT. min_stomate)THEN

             ! Count the number of errors
             err = err + 1

             ! Give extra details of where the mass balance
             ! problem occurs. Stop the model if the user
             ! asked for it
             IF(err_act.GT.1) THEN
                WRITE(numout,*) 'Checking mass balance of a single ipts x ivm'
                WRITE(numout,*) 'Error: mass balance is not closed in ', &
                     TRIM(check_point), ' for ', var_name
                WRITE(numout,*) ' Difference is, ', closure_intern(iele)
                WRITE(numout,*) ' pool_end,pool_start: ', &
                     pool_end(iele), pool_start(iele)
                CALL flush(numout)  
                ! Exit or write warning depening on plev
                CALL ipslerr_p(plev, TRIM(check_point), &
                     'Mass balance error','','')
             ENDIF

          ENDIF

!!$       !+++CHECK+++
!!$       ! A lot of output - only write a message if their is a problem
!!$       ! One message per pixel. The err condition is needed because under
!!$       ! err_act = 2 the model will write an error to the out_orchidee file
!!$       ! but will not be stopped. Therefore it should be checked whether we
!!$       ! have an error or not
!!$       IF (err_act.GT.1 .AND. err.EQ.0) THEN
!!$          WRITE(numout,*) 'Mass balance closure in ', &
!!$               TRIM(check_point), ' for ',var_name
!!$       ENDIF
!!$       !++++++++++

       ENDDO ! nelements


 END SUBROUTINE check_mass_balance_point


!! ================================================================================================================================
!! FUNCTION     : intermediate_mass_balance_check
!!
!>\BRIEF        
!!
!! DESCRIPTION : check whether the mass balance is closed in "the middle" of a stomate_growth_fun_all.
!! Might be difficult to use the exact same code in other subroutines 
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

 SUBROUTINE intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
               circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
               resp_maint, resp_growth, check_intern_init, ipts, j, label, test_type)
   
   !! 0.1 Input variables 
   INTEGER(i_std), INTENT(in)                        :: npts               !! Domain size (unitless)
   REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)     :: circ_class_biomass !! Biomass components of the model tree  
                                                                           !! within a circumference class
                                                                           !! class @tex $(g C ind^{-1})$ @endtex
   REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: circ_class_n       !! Number of individuals in each circ class
                                                                           !! @tex $(m^{-2})$ @endtex
   REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: pool_start         !! Start of this routine 
                                                                           !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
   REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max          !! PFT "Maximal" coverage fraction of a PFT 
                                                                           !! @tex $(m^2 m^{-2})$ @endtex
   REAL(r_std), DIMENSION(:,:), INTENT(in)           :: bm_alloc_tot       !! Allocatable biomass for the whole plant
                                                                           !! @tex $(gC.m^{-2})$ @endtex
   REAL(r_std), DIMENSION(:,:), INTENT(in)           :: gpp_daily          !! PFT gross primary productivity 
                                                                           !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
   REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: atm_to_bm          !! Nitrogen and carbon which is added to the ecosystem to 
                                                                           !! support vegetation growth (gC or gN/m2/day)
   REAL(r_std), INTENT(in)                           :: dt                 !! Time step of the simulations for stomate (days)
   REAL(r_std), DIMENSION(:,:), INTENT(in)           :: resp_maint         !! PFT maintenance respiration 
                                                                           !! @tex $(gC.m^{-2}dt^{-1})$ @endtex    
   REAL(r_std), DIMENSION(:,:), INTENT(in)           :: resp_growth        !! PFT growth respiration 
                                                                           !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
   REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)       :: check_intern_init  !! Contains the components of the internal
                                                                           !! mass balance chech for this routine
                                                                           !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
   INTEGER(i_std), INTENT(in)                        :: ipts, j            !! indices
   CHARACTER(len=*), INTENT(in)                      :: label              !! label to identify the mass balance check

   !! 0.2 Output variables
   
   !! 0.3 Modified variables
   REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: pool_end           !! Start and end pool of this routine 
                                                                           !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
   CHARACTER(*), INTENT(in)                          :: test_type          !! Which test do we want to run on veget_max
   
   !! 0.4 Local variables
   INTEGER(i_std)                                    :: ipar, iele, icir   !! Indices
   INTEGER(i_std)                                    :: imbc,ivm           !! Indices
   CHARACTER(len=80)                                 :: message            !! Label for error message
   REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements):: check_intern       !! Contains the components of the internal
                                                                           !! mass balance chech for this routine
                                                                           !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
   REAL(r_std), DIMENSION(npts,nvm,nelements)        :: closure_intern     !! Check closure of internal mass balance
                                                                           !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
   
!_ ================================================================================================================================

   SELECT CASE(test_type)
   CASE('pft')
      ! Calculate pool_end
      DO ipar = 1,nparts
         DO iele = 1,nelements
            DO icir = 1,ncirc
               pool_end(:,:,iele) = pool_end(:,:,iele) + &
                    (circ_class_biomass(:,:,icir,ipar,iele) * &
                    circ_class_n(:,:,icir) * veget_max(:,:))
            ENDDO
         ENDDO
      ENDDO
      
      ! Calculate mass balance 
      ! Specific processes for carbon
      check_intern(:,:,:,:) = zero
      check_intern(:,:,iatm2land,icarbon) = &
           check_intern_init(:,:,iatm2land,icarbon) + &
           (gpp_daily(:,:) + atm_to_bm(:,:,icarbon)) * &
           dt * veget_max(:,:)
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
            closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                 check_intern(:,:,imbc,iele)
         ENDDO
      ENDDO
      
      ! The call to the intermediate mass balance check
      ! is located outside an nvm-loop. All PFTs should be
      ! have calculated. Check all PFTs. In this call j is
      ! not used. Its value doesn matter.
      message = 'Intermediate mass balance check '//TRIM(label)
      WRITE(numout,*) message, ' for all pfts'
      CALL check_mass_balance(message, closure_intern, npts, pool_end, &
           pool_start, veget_max, test_type, j)

   CASE ('ipft')
      
      ! Calculate pool_end
      DO ipar = 1,nparts
         DO iele = 1,nelements
            DO icir = 1,ncirc
               pool_end(ipts,j,iele) = pool_end(ipts,j,iele) + &
                    (circ_class_biomass(ipts,j,icir,ipar,iele) * &
                    circ_class_n(ipts,j,icir) * veget_max(ipts,j))
            ENDDO
         ENDDO
      ENDDO
      
      ! Calculate mass balance 
      ! Specific processes for carbon
      check_intern(:,:,:,:) = zero
      check_intern(ipts,j,iatm2land,icarbon) = &
           check_intern_init(ipts,j,iatm2land,icarbon) + &
           (gpp_daily(ipts,j) + atm_to_bm(ipts,j,icarbon)) * &
           dt * veget_max(ipts,j)
      check_intern(ipts,j,iatm2land,initrogen) = &
           check_intern_init(ipts,j,iatm2land,initrogen) + & 
           atm_to_bm(ipts,j,initrogen) * dt * veget_max(ipts,j)
      check_intern(ipts,j,iland2atm,icarbon) = -un * &
           (resp_maint(ipts,j) + resp_growth(ipts,j)) * &
           veget_max(ipts,j)
      
      ! Common processes for icarbon and initrogen
      DO iele=1,nelements
         check_intern(ipts,j,ipoolchange,iele) = -un * (pool_end(ipts,j,iele) - &
              pool_start(ipts,j,iele))
      ENDDO
   
      closure_intern = zero
      DO imbc = 1,nmbcomp
         DO iele=1,nelements
            closure_intern(ipts,j,iele) = closure_intern(ipts,j,iele) + &
                 check_intern(ipts,j,imbc,iele)
         ENDDO
      ENDDO

      ! the call is within a PFT-loop. Check for one specific PFT
      message = 'Intermediate mass balance check '//TRIM(label)
      WRITE(numout,*) message, ' for pft, ', j
      CALL check_mass_balance(message, closure_intern, npts, pool_end, &
           pool_start, veget_max, test_type, j)

   CASE DEFAULT

      CALL flush(numout)
      CALL ipslerr_p(3,'intermediate mass balance check', &
           'test_type not correcty defined','case not specified','')

   END SELECT

 END SUBROUTINE intermediate_mass_balance_check 

 
!! ================================================================================================================================
!! FUNCTION     : check_area_change
!!
!>\BRIEF        
!!
!! DESCRIPTION : check whether the change in vegetation is cancelled out by the 
!!               change in non-biological land covers
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
 
   SUBROUTINE check_area_change(check_point, npts, frac_nobio, frac_nobio_new, loss_gain)

    !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: npts                  !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: frac_nobio            !! Fraction of grid cell covered by lakes, land 
                                                                             !! ice, cities, ... before land cover change (unitless, 0-1)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: frac_nobio_new        !! Fraction of grid cell covered by lakes, land 
                                                                             !! ice, cities, ... after land cover change (unitless, 0-1) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: loss_gain             !! losses and gains due to LCC distributed over all
                                                                             !! age classes and thus taking the age-classes into 
                                                                             !! account (unitless, 0-1)
    CHARACTER(*),INTENT(in)                         :: check_point           !! A flag to indicate at which
                                                                             !! point in the code we're doing
                                                                             !! this check

  !! 0.2 Output variables

  !! 0.3 Modified variables

  !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts                   !! Index
    INTEGER(i_std), DIMENSION(npts)                  :: count                  !! counter
    REAL(r_std)                                      :: total                  !! temporary variable for total area of a pixel (unitless, 0-1)
    REAL(r_std), DIMENSION(npts)                     :: change_nobio           !! Change in the non biological fraction in a pixel (0-1, unitless)

!_ ================================================================================================================================

    ! Quality check. It is now expected that the changes cancel out each other
    ! and are exactley zero.
    count(:) = zero
    change_nobio(:) = SUM(frac_nobio_new(:,:) - frac_nobio(:,:),2)
    WHERE(ABS(SUM(loss_gain(:,:),2)+change_nobio(:)).GT.10*EPSILON(un))
       count(:) = un
    END WHERE

    ! Error checking
    IF(SUM(count(:)).GT.zero)THEN

       ! Refine the error message
       DO ipts = 1,npts
          IF (ABS(SUM(loss_gain(ipts,:))+change_nobio(ipts)).GT.10*EPSILON(un)) THEN
             total = ABS(SUM(loss_gain(ipts,:))+change_nobio(ipts))
             WRITE(numout,*) 'ipts, change_nobio, loss_gain, mismatch, ', &
                  ipts, change_nobio(ipts), SUM(loss_gain(ipts,:)), total
             WRITE(numout,*) 'ipts, frac_nobio_new, frac_nobio. ', &
                  ipts, SUM(frac_nobio_new(ipts,:)), SUM(frac_nobio(ipts,:))
          END IF
       END DO
       WRITE(numout,*) 'Found ,',SUM(count(:)), &
            ' pixels for which the total change differs from zero (=0)'
       CALL ipslerr_p(3,TRIM(check_point), &
            'The change in surface areas do not cancel out to zero', '','')

    END IF

  END SUBROUTINE check_area_change


!! ================================================================================================================================
!! FUNCTION     : check_pixel_area
!!
!>\BRIEF        
!!
!! DESCRIPTION : check whether the vegetation and non-biological fractions
!!               still add up to one (=1).
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE check_pixel_area(check_point, npts, veget_max, frac_nobio)

  !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: npts                  !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: veget_max             !! Cover fraction of the vegetation (unitless, 0-1)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: frac_nobio            !! Fraction of grid cell covered by lakes, land 
                                                                             !! ice, cities, ... (unitless, 0-1) 
    CHARACTER(*),INTENT(in)                         :: check_point           !! A flag to indicate at which
                                                                             !! point in the code we're doing
                                                                             !! this check

  !! 0.2 Output variables

  !! 0.3 Modified variables

  !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts                   !! Index
    INTEGER(i_std), DIMENSION(npts)                  :: count                  !! counter

!_ ================================================================================================================================

    ! Quality check. It is still expected that the different vegetation 
    ! fractions in each pixel sums up to exactly one.
    count(:) = zero
    WHERE(ABS(un-SUM(veget_max(:,:),2)-SUM(frac_nobio(:,:),2)).GT.10*EPSILON(un))
       count(:) = un
    END WHERE
    
    ! Error checking
    IF(SUM(count(:)).GT.zero)THEN
       
       ! Refine the error message
       DO ipts = 1,npts
          IF (ABS(un-SUM(veget_max(ipts,:))-SUM(frac_nobio(ipts,:))).GT.10*EPSILON(un)) THEN
             WRITE(numout,*) 'ipts, frac_nobio, veget_max_new, mismatch, ', ipts, &
                  SUM(frac_nobio(ipts,:)), SUM(veget_max(ipts,:)), &
                  un-SUM(veget_max(ipts,:))-SUM(frac_nobio(ipts,:))
          END IF
       END DO     
       WRITE(numout,*) 'Found ,',SUM(count(:)), ' pixels for which the total surface area &
            is not one'
       CALL ipslerr_p(3,TRIM(check_point), &
            'The sum of the vegetation and non-biological fractions',&
            'differs from one','')

    END IF

  END SUBROUTINE check_pixel_area


!! ================================================================================================================================
!! FUNCTION     : check_variable_change
!!
!>\BRIEF        
!!
!! DESCRIPTION : check whether the relationship between veget_max, veget_max_new and 
!!               loss_gain holds
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE check_variable_change(check_point, npts, veget_max, veget_max_new, loss_gain)

  !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: npts                  !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: veget_max             !! Cover fraction of the vegetation before land cover change (unitless, 0-1)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: veget_max_new         !! Cover fraction of the vegetation after land cover change (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: loss_gain             !! losses and gains due to LCC distributed over all
                                                                             !! age classes and thus taking the age-classes into 
                                                                             !! account (unitless, 0-1)
    CHARACTER(*),INTENT(in)                         :: check_point           !! A flag to indicate at which
                                                                             !! point in the code we're doing
                                                                             !! this check

  !! 0.2 Output variables

  !! 0.3 Modified variables

  !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts                   !! Index
    INTEGER(i_std), DIMENSION(npts)                  :: count                  !! counter
    REAL(r_std)                                      :: total                  !! temporary variable for total area of a pixel (unitless, 0-1)
!_ ================================================================================================================================

    ! Quality check. In slowproc veget_max_new was calculated as veget_max
    ! plus loss_gain. If all went well veget_max, loss_gain and
    ! veget_max_new are not changed in between slowproc and
    ! sapiens_lcchange. One of the reason this assumption is valid is
    ! because LCC is calculated at the first day of the year whereas
    ! the other processes that can change veget_max are calculated
    ! the last day of the year. If the above is not true we have to
    ! carefully check what is going on withe veget_max, loss_gain
    ! and/ot veget_max_new.
    count(:) = zero
    WHERE(ABS(SUM(veget_max_new(:,:)-loss_gain(:,:)-veget_max(:,:),2)).GT.10*EPSILON(un))
       count(:) = un
    END WHERE
    
    ! Error checking
    IF(SUM(count(:)).GT.zero)THEN
          
       ! Refine the error message
       DO ipts = 1,npts
          IF (ABS(SUM(veget_max_new(ipts,:)-loss_gain(ipts,:) - &
               veget_max(ipts,:))).GT.10*EPSILON(un)) THEN
             total = ABS(SUM(veget_max_new(ipts,:)-loss_gain(ipts,:) - &
                  veget_max(ipts,:)))
             WRITE(numout,*) 'ipts, veget_max_new, veget_max, loss_gain, mismatch, ', &
                  ipts, SUM(veget_max_new(ipts,:)), SUM(veget_max(ipts,:)), &
                  SUM(loss_gain(ipts,:)), total
          END IF
       END DO     
       WRITE(numout,*) 'Found ,',SUM(count(:)), ' pixels for which the total surface area &
            is not one'
       CALL ipslerr_p(3,TRIM(check_point),&
            'Too big changes in surface area when updating veget_max_new','','')

    END IF
    
  END SUBROUTINE check_variable_change


!! ================================================================================================================================
!! FUNCTION     : check_variable_swap
!!
!>\BRIEF        
!!
!! DESCRIPTION : at one point in sapiens_lcchange veget_max_new should be
!!               equal to veget_max. Check whether this is the case
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE check_variable_swap(check_point, npts, veget_max, veget_max_new)

  !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: npts                  !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: veget_max             !! Cover fraction of the vegetation before land cover change (unitless, 0-1)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: veget_max_new         !! Cover fraction of the vegetation after land cover change (unitless, 0-1)
    CHARACTER(*),INTENT(in)                         :: check_point           !! A flag to indicate at which
                                                                             !! point in the code we're doing
                                                                             !! this check

  !! 0.2 Output variables

  !! 0.3 Modified variables

  !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts                   !! Index
    INTEGER(i_std), DIMENSION(npts)                  :: count                  !! counter
    REAL(r_std)                                      :: total                  !! temporary variable for total area of a pixel (unitless, 0-1)
!_ ================================================================================================================================

    count(:) = zero
    WHERE(ABS(SUM(veget_max(:,:) - veget_max_new(:,:),2)).GT.10*EPSILON(un))
       count(:) = un
    END WHERE

    ! Error checking
    IF(SUM(count(:)).GT.zero)THEN
       
       ! Refine the error message
       DO ipts = 1,npts
          IF (ABS(SUM(veget_max(ipts,:) - veget_max_new(ipts,:))).GT.10*EPSILON(un)) THEN
             total = ABS(SUM(veget_max(ipts,:) - veget_max_new(ipts,:)))
             WRITE(numout,*) 'ipts, veget_max, veget_max_new, mismatch, ', &
                  ipts, SUM(veget_max(ipts,:)), SUM(veget_max_new(ipts,:)), total
          END IF
       END DO
       WRITE(numout,*) 'Found ,',SUM(count(:)), &
            ' pixels for which the updated veget_max differs from veget_max_new'
       CALL ipslerr_p(3,TRIM(check_point),&
            'Updated veget_max but the update differs from veget_max_new', &
            'improve slowproc_adjust_delta_vegetmax','')

    END IF

  END SUBROUTINE check_variable_swap


! ================================================================================================================================
!! SUBROUTINE   : check_biomass_change
!!
!>\BRIEF        Check mass conservation for the biomass pool in sapiens_lcchange
!!
!! DESCRIPTION  : Within sapiens_lcchange, biomass of pfts may change but in a
!!      predictive way depending on whether there was a loss or gain in veget_max.  
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  None
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE check_biomass_change(npts, dt_days, area, circ_class_biomass, &
       circ_class_n, veget_max, loss_gain, atm_to_bm_old, atm_to_bm, &
       harvest_pool, harvest_pool_old, fresh_litter, bm_to_litter_old, &
       biomass_start, biomass_end, n_call)

    !! 0. Variable declaration

    !! 0.1 Input variables
    INTEGER, INTENT(in)                              :: npts                   !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)    :: circ_class_biomass     !! Biomass components of the model tree within a 
                                                                               !! circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: circ_class_n           !! Number of individuals in each circ class
     !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: veget_max              !! "maximal" coverage fraction of a PFT on the ground
                                                                               !! May sum to
                                                                               !! less than unity if the pixel has
                                                                               !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: loss_gain              !! losses and gains due to LCC distributed over all
                                                                               !! age classes and thus taking the age-classes into 
                                                                               !! account (unitless, 0-1)
    INTEGER(i_std), INTENT(in)                       :: n_call                 !! Number indicating whether it is the first or second call
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: atm_to_bm_old          !! Temporary variable to store atm_to_bm at the start of
                                                                               !! the routine @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: atm_to_bm              !! C and N taken from the atmosphere to get C and N to  
                                                                               !! create the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), INTENT(in)                          :: dt_days                !! Time step of vegetation dynamics for stomate (days)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)    :: harvest_pool           !! The wood and biomass that have been
                                                                               !! havested by humans @tex $(gC)$ @endtex 
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)       :: fresh_litter           !! Pool of fresh litter that is being released during @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)            :: area                   !! surface area of the pixel @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)    :: harvest_pool_old       !! Temporary variable to check mass conservation 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)      :: bm_to_litter_old       !! Temporary transfer of biomass to litter
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex


    !! 0.2 Output variables
     
    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)      :: biomass_start          !! biomass before LCC in gC or gN.
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)      :: biomass_end            !! biomass after LCC in gC or gN.

    !! 0.4 Local variables
    INTEGER(i_std)                                   :: iele,icir,ivm,ipts     !! Indices
    REAL(r_std), DIMENSION(nelements)                :: check                  !! Som of pools and fluxes gC or gN pixel-1. 
!_ ================================================================================================================================

     IF (n_call.EQ.1) THEN
        ! Calculate the total biomass at the start of the routine
        DO ipts = 1,npts   
           DO ivm = 1,nvm
              DO iele = 1, nelements
                 DO icir = 1,ncirc
                    ! Multiply with area to make the units easier to understand.
                    ! The unit of biomass_start is gC pixel-1
                    biomass_start(ipts,ivm,iele) = biomass_start(ipts,ivm,iele) + &
                         SUM(circ_class_biomass(ipts,ivm,icir,:,iele)) * &
                         circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm) * &
                         area(ipts) * contfrac(ipts)
                 END DO
              END DO
           END DO
        END DO 
     END IF ! n_call.eq.1

     IF (n_call.EQ.2) THEN
        ! First calculate the total biomass at the end of the routine
        DO ipts = 1,npts   
           DO ivm = 1,nvm
              DO iele = 1, nelements
                 DO icir = 1,ncirc
                    ! Multiply with area to make the units easier to understand.
                    ! The unit of biomass_start is gC pixel-1
                    biomass_end(ipts,ivm,iele) = biomass_end(ipts,ivm,iele) + &
                         SUM(circ_class_biomass(ipts,ivm,icir,:,iele)) * &
                         circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm) * &
                         area(ipts) * contfrac(ipts)
                 END DO
              END DO
 
              ! Check for mass conservation
              check(:) = zero
              IF (loss_gain(ipts,ivm).GE.zero) THEN
                    
                 ! Gain in surface area. Biomass_end includes the newly established
                 ! vegetation. The C and N in this vegetation was taken from the air
                 ! and stored in atm_to_bm. Subtracting atm_to_bm from biomass_end should
                 ! give biomass_start. biopmass_start and biomass_end are in
                 ! gG or gN pixel-1. Atm_to_bm is in gC m-2 dt-1.
                 check(:) = (biomass_start(ipts,ivm,:) - biomass_end(ipts,ivm,:) + &
                      (atm_to_bm(ipts,ivm,:) - atm_to_bm_old(ipts,ivm,:)) * &
                      veget_max(ipts,ivm) * area(ipts) * contfrac(ipts) * dt_days) / &
                      (area(ipts) * contfrac(ipts))

              ELSE
                 
                 ! Loss in surface area. When surface area is lost, the biomass was
                 ! moved into the harvest_pool and fresh_litter. If all these pools are
                 ! accounted for in biomass_end, it should be possible to reconstruct
                 ! biomass_start. biomass_end, biomass_start ar in gC or gN pixel-1.
                 ! fresh_litter is in gC or gN m-2. harvest_pool is in gC or gN pixel-1.
                 ! Note that in case of a species change there will already be some
                 ! biomass in the harvest_pool. In case of a classic land cover change
                 ! this is extremely unlikley for the moment (but this could be the
                 ! case if salavage logging following wind, fires or pests happens in 
                 ! the same time step). Fresh_litter contains part of biomass that was 
                 ! removed as well as bm_to_litter. To check for mass balance conservation
                 ! of the biomass we need to account for the biomass that went into the 
                 ! fresh_litter but not for the bm_to_litter. The latter comes from a
                 ! previous mortality event.
                 check(:) = (biomass_start(ipts,ivm,:) - biomass_end(ipts,ivm,:) - &
                      SUM(fresh_litter(ipts,ivm,:,:)-bm_to_litter_old(ipts,ivm,:,:),1) * &
                      ABS(loss_gain(ipts,ivm)) * area(ipts) * contfrac(ipts) - &
                      SUM(SUM(harvest_pool(ipts,ivm,:,:,:) - harvest_pool_old(ipts,ivm,:,:,:),3),1)) / &
                      (area(ipts) * contfrac(ipts))
                 
              END IF

              ! Check for errors
              DO iele = 1,nelements
                 IF (ABS(check(iele)).GT.min_stomate) THEN
                    WRITE(numout,*) 'ipts, ivm, iele, ', ipts, ivm, iele
                    WRITE(numout,*) 'loss_gain, check, ', loss_gain(ipts,ivm), check(iele)
                    WRITE(numout,*) 'biomass_start, ', biomass_start(ipts,ivm,iele)
                    WRITE(numout,*) 'biomass_end, ', biomass_end(ipts,ivm,iele)
                    WRITE(numout,*) 'diff, ', biomass_end(ipts,ivm,iele)-biomass_start(ipts,ivm,iele)
                    WRITE(numout,*) 'harvest_pool, ', SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2))
                    WRITE(numout,*) 'fresh_litter from biomass, ', &
                         SUM(fresh_litter(ipts,ivm,:,iele)-bm_to_litter_old(ipts,ivm,:,iele),1) * &
                         ABS(loss_gain(ipts,ivm)) * area(ipts) * contfrac(ipts)
                    WRITE(numout,*) 'atm_to_bm, ', (atm_to_bm(ipts,ivm,iele) - atm_to_bm_old(ipts,ivm,iele)) * &
                         veget_max(ipts,ivm) * area(ipts) * contfrac(ipts) * dt_days
                    CALL ipslerr_p(plev,'sapiens_lcchange','Mass balance error for circ_class_biomass,',&
                         'harvest_pool, fresh_litter and/or atm_to_bm','')
                  ENDIF
              END DO
                 
           END DO

        END DO
        
     ENDIF ! n_call.eq.2

   END SUBROUTINE check_biomass_change


!! ================================================================================================================================
!! FUNCTION     : sort_circ_class_biomass
!!
!>\BRIEF        
!!
!! DESCRIPTION : sort the diameters in circ_class_biomass from low to high and
!!               change the order in circ_class_n accordingly
!! 
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : ::circ_class_biomass (gC tree-1),circ_class_n (tree m-2)
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE sort_circ_class_biomass(circ_class_biomass_new,circ_class_n_new)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: circ_class_n_new
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: circ_class_biomass_new

    !! 0.4 Local variables
    INTEGER(i_std)                                     :: inc,icir,jcir        !! Indices
    REAL(r_std), DIMENSION(nparts,nelements)           :: v1                   !! temporay variable for sorting circ_class_biomass
    REAL(r_std)                                        :: v2                   !! temporay variables for sorting circ_class_n
    REAL(r_std)                                        :: sort                 !! flag triggering sorting routine (0 or 1)
    REAL(r_std), DIMENSION(ncirc)                      :: tree_size            !! temporay variables for checking the biomass
                                                                               !! for each circ class @tex ($gC tree^{-1}$) @endtex
    
!_ ================================================================================================================================

    ! Calculate the total biomass in each circumference class
    tree_size(:)=zero
    DO icir=1,ncirc

       ! Only check the order of the carbon pools
       tree_size(icir)=SUM(circ_class_biomass_new(icir,:,icarbon))

    ENDDO

    ! Check whether the order was preserved 
    sort = zero
    DO icir=2,ncirc

       IF(tree_size(icir) .LT. tree_size(icir-1))THEN

          ! The order is not preserved
          sort = un
          EXIT

       ENDIF

    ENDDO

    ! Sort the biomasses. Note that at the same time we also 
    ! have to sort circ_class_n. We made use of the Shell sort
    ! method
    IF (sort==1) THEN

       ! Based on nummerical recipies in Fortran
       inc = 1

       DO
          inc=3*inc+1
          IF (inc .GT. ncirc) EXIT

       ENDDO

       DO
          inc = inc/3

          DO icir=inc+1,ncirc

             v1(:,:) = circ_class_biomass_new(icir,:,:)
             v2 = circ_class_n_new(icir)
             jcir=icir

             DO

                ! Use the carbon pool to sort the circumference classes
                ! but when moving circ classes around, also move the
                ! other elements, i.e. nitrogen.
                IF (SUM(circ_class_biomass_new(jcir-inc,:,icarbon)) .LE. &
                     SUM(v1(:,icarbon))) EXIT
                circ_class_biomass_new(jcir,:,:) = &
                     circ_class_biomass_new(jcir-inc,:,:)
                circ_class_n_new(jcir) = circ_class_n_new(jcir-inc) 
                jcir = jcir-inc
                IF (jcir .LE. inc) EXIT

             ENDDO

             circ_class_biomass_new(jcir,:,:)=v1(:,:)
             circ_class_n_new(jcir)=v2

          ENDDO

          IF (inc .LE. 1) EXIT

       ENDDO

!!$       ! Debug
!!$       WRITE(numout,*) 'sort end - carbon', &
!!$            SUM(circ_class_biomass_new(:,:,icarbon),2),&
!!$            circ_class_n_new(:)
!!$       WRITE(numout,*) 'sort end - nitrogen', &
!!$            SUM(circ_class_biomass_new(:,:,initrogen),2)
!!$       !-

    ENDIF

  END SUBROUTINE sort_circ_class_biomass

!! ================================================================================================================================
!! SUBROUTINE   : Arrhenius_1d
!!
!>\BRIEF         Calculate temperature dependency based on Arrhenius
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: val_arrhenius
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

FUNCTION Arrhenius_1d (kjpindex,temp,ref_temp,energy_act) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
!_ ================================================================================================================================
    
    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))
  END FUNCTION Arrhenius_1d

!!
!================================================================================================================================
!! SUBROUTINE   : Arrhenius_0d
!!
!>\BRIEF         Calculate temperature dependency based on Arrhenius for a
!!               single pixel
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: val_arrhenius
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!_
!================================================================================================================================

FUNCTION Arrhenius_0d (temp,ref_temp,energy_act) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    REAL(r_std),INTENT(in)                        :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)

    !! 0.2 Result

    REAL(r_std)                                   :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
!_
!================================================================================================================================

    val_arrhenius=EXP(((temp-ref_temp)*energy_act)/(ref_temp*RR*(temp)))
  END FUNCTION Arrhenius_0d


!! ================================================================================================================================
!! SUBROUTINE   : Arrhenius_modified_1d
!!
!>\BRIEF         Calculate temperature dependency based on Arrhenius
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: val_arrhenius
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  FUNCTION Arrhenius_modified_1d (kjpindex,temp,ref_temp,energy_act,energy_deact,entropy) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: energy_deact      !! Deactivation Energy (J mol-1)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: entropy           !! Entropy term (J K-1 mol-1)
        
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
!_ ================================================================================================================================    

    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))  &
         * (1. + EXP( (ref_temp * entropy(:) - energy_deact) / (ref_temp * RR ))) &
         / (1. + EXP( (temp(:) * entropy(:) - energy_deact) / ( RR*temp(:))))
         
  END FUNCTION Arrhenius_modified_1d

!! ================================================================================================================================
!! SUBROUTINE   : Arrhenius_modified_0d
!!
!>\BRIEF         Calculate temperature dependency based on Arrhenius
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: val_arrhenius
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  FUNCTION Arrhenius_modified_0d (kjpindex,temp,ref_temp,energy_act,energy_deact,entropy) RESULT ( val_arrhenius )
    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                     :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)    :: temp              !! Temperature (K)
    REAL(r_std), INTENT(in)                       :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: energy_deact      !! Deactivation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: entropy           !! Entropy term (J K-1 mol-1)
        
    !! 0.2 Result

    REAL(r_std), DIMENSION(kjpindex)              :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
!_ ================================================================================================================================    

    val_arrhenius(:)=EXP(((temp(:)-ref_temp)*energy_act)/(ref_temp*RR*(temp(:))))  &
         * (1. + EXP( (ref_temp * entropy - energy_deact) / (ref_temp * RR ))) &
         / (1. + EXP( (temp(:) * entropy - energy_deact) / ( RR*temp(:))))
         
  END FUNCTION Arrhenius_modified_0d

!!
!================================================================================================================================
!! SUBROUTINE   : Arrhenius_modified_nodim
!!
!>\BRIEF         Calculate temperature dependency based on Arrhenius. 
!!               This function is for a single pixel output.
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: val_arrhenius
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!_
!================================================================================================================================

  FUNCTION Arrhenius_modified_nodim (temp,ref_temp,energy_act,energy_deact,entropy) RESULT (val_arrhenius)
    !! 0.1 Input variables

    REAL(r_std),INTENT(in)                        :: temp              !! Temperature (K)
    REAL(r_std),INTENT(in)                        :: ref_temp          !! Temperature of reference (K)
    REAL(r_std),INTENT(in)                        :: energy_act        !! Activation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: energy_deact      !! Deactivation Energy (J mol-1)
    REAL(r_std),INTENT(in)                        :: entropy           !! Entropy term (J K-1 mol-1)

    !! 0.2 Result

    REAL(r_std)                                   :: val_arrhenius     !! Temperature dependance based on
                                                                       !! a Arrhenius function (-)
!_
!================================================================================================================================    

    val_arrhenius=EXP(((temp-ref_temp)*energy_act)/(ref_temp*RR*(temp))) &
         * (1. + EXP( (ref_temp * entropy - energy_deact) / (ref_temp * RR ))) &
         / (1. + EXP( (temp * entropy - energy_deact) / ( RR*temp)))

  END FUNCTION Arrhenius_modified_nodim


!! ================================================================================================================================
!! FUNCTION   : get_printlev
!!
!>\BRIEF        Open out_orchidee text output file, read global PRINTLEV parmeter and return local PRINTLEV_modname
!!
!! DESCRIPTION  : The first time this function is called the out_orchidee output text file is opend and the parameter 
!!                PRINTLEV is read from run.def file. printlev can be accesed each module in ORCHIDEE which 
!!                makes use of constantes_var module. 
!!
!!                This function also reads the parameter PRINTLEV_modname from run.def file. modname is the 
!!                intent(in) character string to this function. If the variable is set in run.def file, the corresponding 
!!                value is returned. Otherwise the value of printlev is returnd as default. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): The local output level for the module set as intent(in) argument.
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  FUNCTION get_printlev ( modname )

    !! 0.1 Input arguments
    CHARACTER(LEN=*), INTENT(IN) :: modname

    !! 0.2 Returned variable
    INTEGER(i_std)               :: get_printlev

    !! 0.3 Local variables
    LOGICAL, SAVE                :: first_printlev=.TRUE.
!OMP THREADPRIVATE(first_printlev)

    !_ ================================================================================================================================

    !! 1.0  Read the global PRINTLEV from run.def and open output file 
    !       This is only done at first call to this function. 
    IF (first_printlev) THEN
       CALL Set_stdout_file('out_orchidee')
       first_printlev=.FALSE.
    END IF

    ! Set default value as the standard printlev
    get_printlev=printlev
    ! Read optional value from run.def file
    CALL getin_p('PRINTLEV_'//modname, get_printlev)

    IF (get_printlev .EQ. 1) THEN
       IF (is_mpi_root.AND.is_omp_root) THEN
          get_printlev=1
       ELSE
          get_printlev=0
       END IF
    END IF

  END FUNCTION get_printlev
!!
!================================================================================================================================
!! FUNCTION   : fm_polynomial_parameter_fitting
!!
!>\BRIEF        program for parameter fitting involving polynomial curve fitting
!! for various forest management scenarios. It calculates and stores polynomial
!! coefficients for different forest management types and diameter ranges. Here's
!! a breakdown of what the subroutine does:
!! DESCRIPTION  : It declares and allocates memory for two arrays:
!!
!!            RDI(2,4): A two-dimensional array used to store data related to forest
!!            management scenarios.
!!            DIA(4): An array used to store diameter-related data.
!!            It initializes poly_rdi to zero. However, the declaration of poly_rdi is missing
!!            in the provided code.
!!
!!            It enters a loop DO ivm=2,nvm where ivm iterates over a range of values from 2
!!            to nvm.
!!
!!            Inside the loop, it checks if is_tree(ivm) is true, indicating that the current
!!            scenario is related to a tree.
!!
!!            If the condition is true, it proceeds to calculate polynomial coefficients for
!!            different forest management types (e.g., NONE, THIN, UNEVEN, COPPICE) within the
!!            context of the current tree (ivm).
!!
!!            For each forest management type, it sets up the DIA and RDI arrays with
!!            appropriate values.
!!
!!            It calls the polyfit function to fit a polynomial curve to the data in the DIA
!!            and RDI arrays and stores the resulting polynomial coefficients in the poly_rdi
!!            array. The specific polynomial degree used in each fitting operation varies
!!            (e.g., 2 or 3), depending on the scenario.
!!
!!            It writes the fitted polynomial coefficients to the output, providing
!!            information about the upper and lower parameters for each forest management type
!!            within the tree scenario.
!!
!!            After completing the loop, the subroutine deallocates the memory used for the
!!            RDI and DIA arrays.
!!
!!            It's worth noting that some details are missing from the provided code, such as
!!            the declaration and allocation of the poly_rdi array and the definitions of
!!            various variables and functions (e.g., nvm, is_tree, largest_tree_dia,
!!            rdi_start, thin_intensity, etc.). These details are essential for a complete
!!            understanding of the code and its functionality.
!!
!!            Overall, the subroutine seems to be part of a larger program that deals with
!!            forest management and parameter fitting, particularly with polynomial curve
!!            fitting for different scenarios related to trees.
!!
!! MAIN OUTPUT VARIABLE(S): poly_rdi
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
  SUBROUTINE fm_polynomial_parameter_fitting()
    REAL(r_std), DIMENSION(:), allocatable     :: DIA 
    REAL(r_std), DIMENSION(:,:), allocatable   :: RDI
    INTEGER(i_std)                             :: ivm
    REAL(r_std)                                :: dia_init_mean
    ALLOCATE(RDI(2,4))
    ALLOCATE(DIA(4))


    poly_rdi = zero

    DO ivm=2,nvm

      IF (is_tree(ivm)) THEN

        RDI = zero
        DIA = zero

        IF (printlev >= 3) WRITE(numout,*) "For PFT : ", ivm, " calculate RDI param" 
   
        !commun
        DIA(1) = zero
        RDI(1,1) = rdi_start(ivm)
        RDI(1,2) = rdi_start(ivm)
        RDI(2,1) = rdi_start(ivm) * (1-thin_intensity(ivm))
        RDI(2,2) = rdi_start(ivm) * (1-thin_intensity(ivm))

        ! None : unmanaged forest
        DIA(4) = final_dia(ivm, ifm_none)*100.0
        DIA(3) = DIA(4)*r_max
        RDI(1,3) = rdi_max(ivm, ifm_none)
        RDI(1,4) = rdi_dia_ref(ivm)
        RDI(2,3) = RDI(1,3) * (1-thin_intensity(ivm))
        RDI(2,4) = RDI(1,4) * (1-thin_intensity(ivm))
 
        ! fitting polynomial parameters
        poly_rdi(ivm, ifm_none, 1, :) = polyfit(DIA([1,3,4]), RDI(1,[1,3,4]), 2)
        poly_rdi(ivm, ifm_none, 2, :) = polyfit(DIA([1,3,4]), RDI(2,[1,3,4]), 2)
 
        IF (printlev >= 3) WRITE(numout,*) "RDI parameter upper for ifm=NONE : ", &
             poly_rdi(ivm, ifm_none, 1, :)
        IF (printlev >= 3) WRITE(numout,*) "RDI parameter low for ifm=NONE : ", &
             poly_rdi(ivm,ifm_none, 2, :)

        ! Thin : thinnning and clearcut even-aged forest
        DIA(2) = dia_boost_end(ivm)
        DIA(4) = final_dia(ivm, ifm_thin)*100.0
        DIA(3) = DIA(4)*r_max
        RDI(1,3) = rdi_max(ivm, ifm_thin)
        RDI(1,4) = rdi_max(ivm, ifm_thin)
        RDI(2,3) = RDI(1,3) * (1-thin_intensity(ivm))
        RDI(2,4) = RDI(1,4) * (1-thin_intensity(ivm))

        ! fitting polynomial parameters
        poly_rdi(ivm, ifm_thin, 1, :) = polyfit(DIA, RDI(1,:), 3)
        poly_rdi(ivm, ifm_thin, 2, :) = polyfit(DIA, RDI(2,:), 3)
  
        IF (printlev >= 3) WRITE(numout,*) "RDI parameter upper for ifm=THIN : ", &
             poly_rdi(ivm,ifm_thin, 1, :)
        IF (printlev >= 3) WRITE(numout,*) "RDI parameter low for ifm=THIN : ", &
             poly_rdi(ivm,ifm_thin, 2, :)

        ! Uneven : Continuous cover forestry
        DIA(4) = final_dia(ivm, ifm_uneven)*100.0
        DIA(3) = DIA(4)*r_max
        RDI(1,3) = rdi_max(ivm, ifm_uneven)
        RDI(1,4) = rdi_max(ivm, ifm_uneven) * (rdi_dia_ref(ivm)/rdi_max(ivm, ifm_thin))
        RDI(2,3) = RDI(1,3) * (1-thin_intensity(ivm))
        RDI(2,4) = RDI(1,4) * (1-thin_intensity(ivm))
        
        ! fitting polynomial parameters   
        poly_rdi(ivm, ifm_uneven, 1, :) = polyfit(DIA, RDI(1,:), 3)
        poly_rdi(ivm, ifm_uneven, 2, :) = polyfit(DIA, RDI(2,:), 3)

        IF (printlev >= 3) WRITE(numout,*) "RDI parameter upper for ifm=UNEVEN : ", &
             poly_rdi(ivm,ifm_uneven, 1, :)
        IF (printlev >= 3) WRITE(numout,*) "RDI parameter low for ifm=UNEVEN : ", &
             poly_rdi(ivm,ifm_uneven, 2, :)

        !+++CHECK+++
        !code still needs to be developed for coppice and short rotation coppice
        !+++++++++++

     END IF
    END DO
    DEALLOCATE(RDI)
    DEALLOCATE(DIA)

  END SUBROUTINE


!!
!================================================================================================================================
!! FUNCTION   : polyfit
!!
!>\BRIEF     :  performing polynomial curve fitting using a least squares
!method. It fits a polynomial of degree d to the input data points (vx, vy).
!!
!! DESCRIPTION  : It initializes and allocates memory for various arrays and
!!variables.
!!
!!                It prepares the Vandermonde matrix X where each column contains powers of the
!!                input values vx. This matrix is used to represent the polynomial equations.
!!
!!                It initializes the polyfit array to zero. This array will store the coefficients
!!                of the polynomial fit.
!!
!!                It calculates the transpose of the matrix X and the matrix product XTX (X
!!                transpose times X).
!!
!!                It performs LU decomposition on the XTX matrix using the ludcmp subroutine. LU
!!                decomposition is a method for solving linear systems of equations.
!!
!!                It enters a loop where it iteratively improves the polynomial coefficients using
!!                the mprove subroutine and calculates the difference poly_diff between the
!!                current polynomial fit and the data points. The loop continues until poly_diff
!!                is below a specified threshold (res) or until a maximum number of iterations is
!!                reached.
!!
!!               If the loop exits without convergence (i.e., the maximum number of iterations is
!!               reached), it writes a message indicating that there was no convergence.
!!
!!               Finally, it deallocates the dynamically allocated arrays and returns the
!!               polynomial coefficients in the polyfit array.
!!
!!               The code also mentions a gauss_jordan_method subroutine, which is commented out
!!               and not used in the code.
!!
!!               It's important to note that this code assumes the availability of certain
!!               subroutines (ludcmp, mprove, and test_param)  These subroutines contain 
!!               the specific numerical algorithms for LU decomposition, 
!!               polynomial coefficient improvement, and testing convergence.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): polyfit
!!
!! REFERENCE(S) : ludcmp, mprove, and test_param
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
  FUNCTION polyfit(vx, vy, d)
        integer(i_std), intent(in)               :: d
        real(r_std), dimension(d+1)              :: polyfit
        real(r_std), dimension(:), intent(in)    :: vx, vy
       
        real(r_std), dimension(:,:), allocatable :: X
        real(r_std), dimension(:,:), allocatable :: XT
        real(r_std), dimension(:,:), allocatable :: XTX
       
        integer(i_std)                           :: i, j 
        integer(i_std)                           :: n, iter
        integer, dimension(:), allocatable       :: ipiv
        real(r_std)                              :: dd, res, poly_diff
        

        n = d+1
        res = 1e-5
        poly_diff = 1.0
        allocate(ipiv(n))
        allocate(XT(n, size(vx)))
        allocate(X(size(vx), n))
        allocate(XTX(n, n))
       
        ! prepare the matrix
        do i = 0, d
           do j = 1, size(vx)
              X(j, i+1) = vx(j)**i
           end do
        end do
        polyfit = 0.0
        XT  = transpose(X)
        XTX = matmul(XT, X)

        ! LU decomposition
        call ludcmp(XTX,ipiv,dd)
        iter = 0
        do while (poly_diff > res .AND. iter < 1e6)
          iter = iter + 1
          call mprove(X, XTX, ipiv, vy, polyfit)
          call test_param(polyfit, vx, vy, d, poly_diff)
          if (iter == 1e6) WRITE(numout,*) "no convergance at 1M eteration diff=", poly_diff
        end do

        IF (printlev > 3) WRITE(numout,*) "final : diff=", poly_diff, "at ", iter, "iterations" 
       
        deallocate(ipiv)
        deallocate(X)
        deallocate(XT)
        deallocate(XTX)
      
   END FUNCTION

!!
!================================================================================================================================
!! FUNCTION   : test_parma
!!
!>\BRIEF      : calculate the difference between the predicted values of a
!!              polynomial fit and the actual data points. 
!!
!! DESCRIPTION  : p: An array of coefficients representing the polynomial fit of
!!                   degree d.
!!                x: An array of input data points.
!!                y: An array of actual data points.
!!                d: The degree of the polynomial.
!!                diff: An input/output variable used to store the calculated difference.
!!                It initializes an array vy of the same size as y to store the predicted values.
!!
!!                It loops through each data point in the x array and calculates the predicted
!!                value vy(j) for each data point using the polynomial coefficients p. The
!!                predicted value is computed by evaluating the polynomial at the given input
!!                value x(j).
!!
!!                It calculates the difference diff between the predicted values vy and the actual
!!                data points y using the root mean square error (RMSE) formula:
!!
!!               diff = SQRT(SUM((vy - y)**2) / size(x))
!!               This formula computes the square root of the average of the squared differences
!!               between the predicted and actual values. It's a measure of how well the
!!               polynomial fit matches the data points, with lower values indicating a better
!!               fit.
!!
!!               Overall, this subroutine helps assess the quality of the polynomial fit by
!!               calculating the RMSE between the fitted curve and the actual data. The diff
!!               variable will hold this RMSE value after the subroutine is called. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): diff
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
  SUBROUTINE test_param(p, x, y, d, diff)
        integer(i_std), intent(in)                  :: d
        real(r_std), dimension(:), intent(in)       :: x, y
        real(r_std), dimension(d+1), intent(in)     :: p

       real(r_std), intent(inout)                   :: diff
       real(r_std), dimension(size(y))              :: vy
       integer(i_std)                               :: i,j,n
       
       vy = 0.0
       n=d+1
       do j=1, size(x)
         vy(j) = p(1)
         do i=2, n
           vy(j) = vy(j) + p(i) * x(j)**(i-1) 
         end do
       end do

       diff = SQRT(SUM((vy-y)**2)/size(x))

  END SUBROUTINE

!!
!================================================================================================================================
!! FUNCTION   : ludcmp
!!
!>\BRIEF      : Perform a LU matrix decomposition
!!
!! DESCRIPTION  :  Given an N × N input matrix a, this routine replaces it by the LU
!!                 decomposition of a
!!                 rowwise permutation of itself. On output; indx is an
!!                 output vector of length N that records the row permutation effected by
!!                 the partial pivoting;
!!                 d is output as ±1 depending on whether the number of row interchanges
!!                 was even or odd,
!!                respectively. This routine is used in combination with lubksb to solve
!!                linear equations or invert a matrix.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): a, indx, d
!!
!! REFERENCE(S) : swap, imaxloc, outerprod
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
   SUBROUTINE ludcmp(a,indx,d)
        REAL(r_std), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER(i_std), DIMENSION(:), INTENT(OUT)  :: indx
        REAL(r_std), INTENT(OUT)                   :: d
        REAL(r_std), DIMENSION(size(a,1)) :: vv !vv stores the implicit scaling of each row.
        REAL(r_std), PARAMETER :: TINY=1.0e-20 !A small number.
        INTEGER(i_std) :: j,n,imax
        n=size(a,1)
        d=un !No row interchanges yet.
        vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling
        !if (any(vv == 0.0)) call nrerror(’singular matrix in ludcmp’)
        vv=un/vv !Save the scaling.
        do j=1,n
            imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
            if (j /= imax) then !Do we need to interchange rows?
                call swap(a(imax,:),a(j,:)) !Yes, do so...               
                d=-d !...and change the parity of d.
                vv(imax)=vv(j) !Also interchange the scale factor.
            end if
            indx(j)=imax
            if (a(j,j) == 0.0) a(j,j)=TINY
            !If the pivot element is zero the matrix is singular (at least to
            !the precision of the algorithm). 
            !For some applications on singular matrices, it is desirable to
            !substitute TINY for zero.
            a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
            a(j+1:n,j+1:n)=a(j+1:n,j+1:n) - outerprod(a(j+1:n,j),a(j,j+1:n))
            !Reduce remaining submatrix.
        end do
     END SUBROUTINE ludcmp

!!
!================================================================================================================================
!! FUNCTION   : lubksb
!!
!>\BRIEF        Solves the set of N linear equations A · X = B
!!
!! DESCRIPTION  : Here the N × N matrix a is input, not
!! as the original matrix A, but rather as its LU decomposition,
!! determined by the routine
!! ludcmp. indx is input as the permutation vector of length N returned by
!! ludcmp. b is
!! input as the right-hand-side vector B, also of length N, and returns
!! with the solution vector
!! X. a and indx are not modified by this routine and can be left in place
!! for successive calls
!! with different right-hand sides b. This routine takes into account the
!! possibility that b will
!! begin with many zero elements, so it is efficient for use in matrix inversion.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): b
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
     SUBROUTINE lubksb(a,indx,b)
        REAL(r_std), DIMENSION(:,:), INTENT(IN)  :: a
        INTEGER(i_std), DIMENSION(:), INTENT(IN) :: indx
        REAL(r_std), DIMENSION(:), INTENT(INOUT) :: b
        INTEGER(i_std) :: i,n,ii,ll
        REAL(r_std) :: summ
        n=size(a,1)
        ii=0  !When ii is set to a positive value, it will become the index of
              !the first nonvanishing element of b. We now do
              !the forward substitution. The only new
              !wrinkle is to unscramble the permutation as we go.
        do i=1,n
            ll=indx(i)
            summ=b(ll)
            b(ll)=b(i)
            if (ii /= 0) then
                summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
            else if (summ /= 0.0) then
                ii=i !A nonzero element was encountered, so from now on we will
            end if !have to do the dot product above.
            b(i)=summ
        end do
        do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
            b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
        end do
     END SUBROUTINE lubksb

!!
!================================================================================================================================
!! FUNCTION   : mprove
!!
!>\BRIEF        Improves a solution vector x of the linear set of equations A· X= B.
!!
!! DESCRIPTION  : The N ×N matrix a
!! and the N-dimensional vectors b and x are input. Also input is alud,
!! the LU decomposition
!! of a as returned by ludcmp, and the N-dimensional vector indx also
!! returned by that
!! routine. On output, only x is modified, to an improved set of values.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): x
!!
!! REFERENCE(S) : lubksb
!!
!! FLOWCHART    :
!! \n
!_
!================================================================================================================================
     SUBROUTINE mprove(a,alud,indx,b,x)
        REAL(r_std), DIMENSION(:,:), INTENT(IN)  :: a,alud
        INTEGER(i_std), DIMENSION(:), INTENT(IN) :: indx
        REAL(r_std), DIMENSION(:), INTENT(IN)    :: b
        REAL(r_std), DIMENSION(:), INTENT(INOUT) :: x
        INTEGER(i_std) :: ndum
        REAL(r_std), DIMENSION(size(a,1)) :: r
        ndum=size(a,1)
        r=matmul(real(a,r_std),real(x,r_std))-real(b,r_std)
        !Calculate the right-hand side, accumulating the residual in double
        !precision.
        call lubksb(alud,indx,r) !Solve for the error term,
        x=x-r !and subtract it from the old solution.
     END SUBROUTINE mprove

!!
!================================================================================================================================
!! FUNCTION   : swap
!!
!>\BRIEF        simple routine for swapping the contents of two one-dimensional arrays a and b.
!!
!! \n
!_
!================================================================================================================================
     SUBROUTINE swap(a,b)
        REAL(r_std), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(r_std), DIMENSION(SIZE(a))          :: dum
        dum=a
        a=b
        b=dum
     END SUBROUTINE swap

!!
!================================================================================================================================
!! FUNCTION   : imaxloc
!!
!>\BRIEF       function that returns the index of the maximum value in a one-dimensional array arr
!!
!! \n
!_
!================================================================================================================================
     FUNCTION imaxloc(arr)
        !Index of maxloc on an array.
        REAL(r_std), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(i_std)                        :: imaxloc
        INTEGER(i_std), DIMENSION(1)          :: imax
        imax=maxloc(arr(:))
        imaxloc=imax(1)
     END FUNCTION imaxloc

!!
!================================================================================================================================
!! FUNCTION   : outerprod
!!
!>\BRIEF      function that calculates the outer product of two one-dimensional
!!            arrays a and b to create a two-dimensional array.
!!
!! \n
!_
!================================================================================================================================
     FUNCTION outerprod(a,b)
        REAL(r_std), DIMENSION(:), INTENT(IN)   :: a,b
        REAL(r_std), DIMENSION(size(a),size(b)) :: outerprod
        outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
     END FUNCTION outerprod

END MODULE function_library
