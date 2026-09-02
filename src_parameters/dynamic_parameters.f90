! =================================================================================================================================
! MODULE       : dynamic_parameters
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module initializes dynamic parameters.
!!
!!\n DESCRIPTION:  This module allocates and initializes dynamic parameters in function of the number of pfts, 
!!                 number of pixels, and an initial value of the parameters. \n
!!                 The number of pixels is read in the driver
!!                 The number of PFTs is read in control.f90 (subroutine control_initialize). \n
!!                 Then we can initialize the parameters. \n
!!
!! RECENT CHANGE(S): Sebastiaan Luyssaert 2023: Added dynamic parameters
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_parameters/pft_parameters.f90 $
!! $Date: 2023-09-18 17:37:23 +0200 (Mon, 18 Sep 2023) $
!! $Revision: 8176 $
!! \n
!_ ================================================================================================================================

MODULE dynamic_parameters

  USE dynamic_parameters_var
  USE pft_parameters_var
  USE constantes_mtc
  USE constantes_var
  USE ioipsl
  USE ioipsl_para
  USE defprec
  
  IMPLICIT NONE      

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_main
!!
!>\BRIEF          This subroutine initializes all the dynamic parameters in function of the
!! number of vegetation types and pixel chosen by the user.
!!
!! DESCRIPTION  : 
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

  SUBROUTINE dynamic_parameters_main(kjpindex)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    INTEGER(i_std), INTENT(in)          :: kjpindex                           !! Number of grid cells in the spatial domain (-)   
    

    !_ ================================================================================================================================ 

    IF(l_first_dynamic_parameters) THEN

       !! First time step
       IF(printlev>=3) THEN
          WRITE(numout,*) 'l_first_dynamic_parameters: allocate dynamic parameters'
          CALL flush(numout)
       END IF
          
       !! Memory allocation for the dynamic parameters
       CALL dynamic_parameters_alloc(kjpindex)


       !! Initialize dynamic parameters
       CALL dynamic_parameters_init(kjpindex)

       IF (printlev>=3) THEN
          WRITE(numout,*) 'first_dynamic_parameters done'
          CALL flush(numout)
       END IF
       
       !! Reset flag
       l_first_dynamic_parameters = .FALSE.

    ELSE 

       RETURN

    ENDIF !(l_first_dynamic_parameters)

  END SUBROUTINE dynamic_parameters_main


!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_init 
!!
!>\BRIEF          This subroutine initializes the dynamic parameters by the default values
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S):  
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE dynamic_parameters_init(kjpindex)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)          :: kjpindex                           !! Number of grid cells in the spatial domain (-) 

    !! 0.4 Local variables
    INTEGER(i_std)                      :: jv                                 !! Index (unitless)
    
    !_ ================================================================================================================================ 

    ! 1.2 For sechiba parameters
    IF (ok_sechiba) THEN

       pipe_tune2(:,ibare_sechiba) = zero
       
       ! pipe_tune2 will be calculated as a function
       ! of precipitation. Because it takes a year to make a first
       ! estimate of the annual precipitation, pipe_tune2 needs to
       ! be initialized with a reasonable small value. If not the
       ! model may crash. In cases where there is no rain in on the
       ! first day pipe_tune2 would equal to zero. 
       DO jv = 2,nvm
          pipe_tune2(:,jv) = min_pipe_tune2(jv)
       END DO
       
       IF (ok_stomate) THEN
       
          IF (.NOT. test_dynamic_alpha_self_thin) THEN
             ! When using a static alpha_self_thinning, the value is initialized
             ! here and never adjusted afterwards.
             DO jv = 1,nvm
                alpha_self_thinning(:,jv) = alpha_self_thinning_mtc(pft_to_mtc(jv))
             END DO
          ELSE
             ! When using the dynamic alpha_self_thinning, the value is initialized
             ! here but it will get recalculated in stomate_lpj before it is used.
             ! As such it does not matter how it is initilaized. Both cases were
             ! nevertheless separated and initialized for clarity and consistency.
             DO jv = 1,nvm
                alpha_self_thinning(:,jv) = ref_alpha_self_thin_mtc(pft_to_mtc(jv))
             END DO
          END IF

       END IF ! ok_stomate

    END IF ! ok_sechiba

  END SUBROUTINE dynamic_parameters_init


!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_alloc
!!
!>\BRIEF         This subroutine allocates memory needed for the dynamic parameters 
!! in function of the flags activated.  
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

  SUBROUTINE dynamic_parameters_alloc(kjpindex)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)          :: kjpindex                           !! Number of grid cells in the spatial domain (-)
    
    !! 0.4 Local variables
    LOGICAL                             :: l_error                            !! Diagnostic boolean for error allocation (true/false) 
    INTEGER                             :: ier                                !! Return value for memory allocation (0-N, unitless)

    !_ ================================================================================================================================

    !! Parameters used if ok_sechiba only
    IF ( ok_sechiba ) THEN

       l_error = .FALSE.
       ALLOCATE(pipe_tune2(kjpindex,nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pipe_tune2. We stop. We need kjpindex * nvm words = ',kjpindex, nvm
          CALL ipslerr_p (3,'dynamic_parameters','dynamic_parameters_alloc','','')
       END IF

       ! Should depend on space and time. Not entirely clear whether this should be
       ! condidered a variable or a parameter. Because it can receive fixed values
       ! from constantes_mtc it is initialized as a parameter, but its usage is
       ! closer to that of a variable.
       ALLOCATE(alpha_self_thinning(kjpindex,nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alpha_self_thinning. We stop. We need kjpindex * nvm words = ', kjpindex, nvm
          CALL ipslerr_p (3,'dynamic_parameters','dynamic_parameters_alloc','','')
       END IF

       ! Should depend on space and time. Not entirely clear whether this should be
       ! condidered a variable or a parameter. Because it can receive fixed values
       ! from constantes_mtc it is initialized as a parameter, but its usage is
       ! closer to that of a variable.
       ALLOCATE(pre_indust_ref_gpp(kjpindex,nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pre_indust_ref_gpp. We stop. We need kjpindex * nvm words = ',kjpindex, nvm
          CALL ipslerr_p (3,'dynamic_parameters','dynamic_parameters_alloc','','')
       END IF
       
    END IF

  END SUBROUTINE dynamic_parameters_alloc


!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_config_sechiba
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

  SUBROUTINE dynamic_parameters_config_sechiba()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variable
    INTEGER(i_std)                                      :: jv                 !! Index (unitless)
    REAL(r_std), DIMENSION(nvm)                         :: temp               !! temporary variable enabling to use getin for variables
                                                                              !! with other dimensions than nvm
    !_ ================================================================================================================================ 

    !Config Key   = PIPE_TUNE2 
    !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 40., 40., 40., 40., 40., 40., 40., 40., 0., 0., 0., 0.  
    !Config Help  = 
    !Config Units = [-]
    temp(:) = zero
    CALL getin_p('PIPE_TUNE2',temp)

    ! The dimensions do not match. getin_p was made to put a scalar
    ! into an array(nvm). pipe_tune2 is now an array (npts,nvm).
    ! Use temp as an intermediate
    DO jv = 1,nvm
       IF (temp(jv).GT.zero) THEN
          pipe_tune2(:,jv) = temp(jv)
       END IF
    END DO
    
  END SUBROUTINE dynamic_parameters_config_sechiba

  
!! ================================================================================================================================

  SUBROUTINE dynamic_parameters_config_stomate()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variable
    INTEGER(i_std)                                      :: jv                 !! Index (unitless)
    REAL(r_std), DIMENSION(nvm)                         :: temp               !! temporary variable enabling to use getin for variables
                                                                              !! with other dimensions than nvm
    !_ ================================================================================================================================ 

    !Config Key   = ALPHA_SELF_THINNING
    !Config Desc  = alpha coefficient of the self thinning relationship
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 2827, 3700, 1348, 1220, 2000, 2827, 800, 800, undef, undef, undef, undef
    !Config Help  = Coefficient of the self-thinning relationship D=alpha*N^beta 
    !               estimated from German, French, Spanish and Swedish 
    !               forest inventories
    !Config Units = [-]
    temp(:) = zero
    CALL getin_p("ALPHA_SELF_THINNING",temp)

    ! The dimensions do not match. getin_p was made to put a scalar
    ! into an array(nvm). alpha_self_thinning is now an array (npts,nvm).
    ! Use temp as an intermediate
    DO jv = 1,nvm
       IF (temp(jv).GT.zero) THEN
          alpha_self_thinning(:,jv) = temp(jv)
       END IF
    END DO
    
  END SUBROUTINE dynamic_parameters_config_stomate

!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_clear
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

  SUBROUTINE dynamic_parameters_clear

    l_first_dynamic_parameters = .TRUE.

    IF (ALLOCATED(pipe_tune2)) DEALLOCATE(pipe_tune2)
    IF (ALLOCATED(alpha_self_thinning)) DEALLOCATE(alpha_self_thinning)
    IF (ALLOCATED(pre_indust_ref_gpp)) DEALLOCATE(pre_indust_ref_gpp)

  END SUBROUTINE dynamic_parameters_clear


!! ================================================================================================================================
!! SUBROUTINE   : dynamic_parameters_calculate
!!
!>\BRIEF         This subroutine calculates dynamic parameters 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : Kempes CP, West GB, Crowell K, Girvan M (2011) Predicting Maximum Tree Heights
!!                and Other Traits from Allometric Scaling and Resource Limitations. PLoS ONE
!!                6(6): e20551. doi:10.1371/journal.pone.0020551
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE dynamic_parameters_calculate(kjpindex, precip_longterm)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex                  !! Number of grid cells in the spatial domain (-) 
    REAL(r_std), DIMENSION(:), INTENT(in)               :: precip_longterm           !! longterm annual precipitation sum 
                                                                                     !! @tex ($mm year^{-1}$) @endtex

    !! 0.4 Local variables
    INTEGER(i_std)                                      :: jv, ipts                  !! Index (unitless)
    !_ ================================================================================================================================

    ! The parameter still has the nvm dimension but that is just to allow for
    ! both approaches. If the dynamic approach is kept, the pft-dimension could be
    ! removed. Doing so would imply that at a given pixel, all PFTs can reach the
    ! same tree height. This is an acceptable large scale assumption as demonstrated
    ! by Kempes et al 2011 for the USA.

    ! The model itself is taken from kempes et al 2011 which states that transpiration
    ! of an individual tree (Q0) scales to tree height:
    ! Q0=beta2*h^neta2 (Eq. 3 in Kempes)
    ! The units given in the paper for Eq.3 make sense. The root radius also scales
    ! to tree height:
    ! rroot=beta3**(1/4)*h (Eq. 4 in Kempes)
    ! No units for beta 3 given in the paper -> rroot is thus in cm.
    ! The water available to an individual tree (Qp) can be calculated as:
    ! Qp=gamma*pi*(rroot**2)*(rain/1000)
    ! It is written that rain is in m/year. The paper does not give
    ! units for the other variables. gamma may convert everything but
    ! it should then have units. It doesn't. Whithout unit conversion
    ! the calculations don't make any sense. The first author was contacted (in 2023)
    ! to clarify the units of the parameters but did not reply to my [SL] email.
    ! For this reason it was decided to fit the parameters of this model ourself.
    
    ! The maximum tree height will be reached when Q0 = Qp. Following substitution,
    ! the model can be simplified to:
    ! height = alpha*precip^beta
    ! This is the model that was parameterized using post-processed GEDI for tree
    ! height and CRUJRA for precipitation. The script used can be found at
    ! svn/TOOLS/FIT_TREEHEIGHT_PRECIPITATION/dynamic_tree_height.py
    
    ! The parameter values were obtrained by fitting CRUJRA to GEDI for the year 2000.
    ! When fitting the parameters against the observations we are fitting a relationship
    ! between the 75% percentile of tree height vs precipitation. We don know the
    ! diameter of this 75% percentile of tree height. This is an issue because ORCHIDEE
    ! use diameter-height relationship and in this relationship pipe_tune2 is the height
    ! for a tree with a diameter of 1 m. It was assumed that the 75% percentile height
    ! was observed for trees with a diameter of 50 cm. Hence the fitted height had to
    ! be converted from 0.5m to 1.0m. The value obtained for alpha_height_precip for
    ! the assumed 0.5m was 2.33. Following conversion the value for alpha_height_precip
    ! used in this approach is 3.22. Given the need of this conversion and the underlying
    ! fitting, these values could be optimized. This conversion is documented in
    ! svn/TOOLS/FIT_TREEHEIGHT_PRECIPITATION/dynamic_tree_height.py 
    
    ! Note that the parameters were fit for annual preciptation. The ORCHIDEE driver
    ! gives precipitation in mm/day. So a conversion is need. Also, the first day,
    ! when vegetation is prescribed, precip_longterm could still be zero. This would
    ! make the model crash, hence, the use of a mimimum value (min_pipe_tune2).
    DO ipts = 1, kjpindex
       DO jv = 2,nvm
          pipe_tune2(ipts,jv) =  MAX(min_pipe_tune2(jv),&
               alpha_height_precip*((precip_longterm(ipts)*365)**beta_height_precip))
       END DO
    END DO

  END SUBROUTINE dynamic_parameters_calculate  

END MODULE dynamic_parameters
