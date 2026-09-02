! =================================================================================================================================
! MODULE       : stomate_stand_structure
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_stand_structure.f90 $
!! $Date: 2023-11-22 20:28:54 +0100 (mer. 22 nov. 2023) $
!! $Revision: 8316 $
!! \n
!_ ================================================================================================================================

MODULE stomate_stand_structure

  ! modules used:

  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE dynamic_parameters
  USE constantes
  USE function_library, ONLY:wood_to_height, wood_to_dia, &
                             wood_to_cv, wood_to_cn_dia

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC stand_structure_clear, derive_biomass_quantities

  ! first call
  LOGICAL, SAVE                                    :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall) 
CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : stand_structure_clear
!!
!>\BRIEF        : Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE stand_structure_clear
    firstcall=.TRUE.
  END SUBROUTINE stand_structure_clear


!! ================================================================================================================================
!! SUBROUTINE 	:derive_biomass_quantities
!!
!>\BRIEF        Use the basal areabiomass and number density to derive various
!!              distributions of the trees in a single grid point for
!!              all vegetation types
!!
!! DESCRIPTION  : I have chosen to do this for a single grid point instead of 
!!                the whole map or a single grid point and single PFT because
!!                of the compromise between speed (subroutine overhead) and
!!                flexibility
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::height_dist, ::diameter_dist, ::cn_area_dist, ::cn_vol_dist
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE derive_biomass_quantities(npts, nvm, ncirc, circ_class_n, &
       circ_class_biomass, values)

  !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    INTEGER,INTENT(IN)                                      :: npts               !! Number of pixels
    INTEGER,INTENT(IN)                                      :: nvm                !! Number of PFT types
    INTEGER,INTENT(IN)                                      :: ncirc              !! Number of circumference classes
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN)           :: circ_class_biomass !! Biomass of the componets of the model  
                                                                                  !! tree within a circumference
                                                                                  !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)               :: circ_class_n       !! Number of trees within each circumference
                                                                                  !! class @tex $(m^{-2})$ @endtex

    !! 0.2 Output variables

    REAL(r_std),DIMENSION(:,:,:,:),INTENT(OUT)              :: values             !! An array which holds data for
                                                                                  !! various canopy parameters 


    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                          :: ipts, ivm, icir, &
                                                                       idist_type !! index (unitless)
!_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) 'Entering derive_biomass_quantities'

    ! zero everything
    values(:,:,:,:)=zero

    DO ipts=1,npts
       DO ivm=1,nvm

          IF (.NOT. is_tree(ivm)) CYCLE

          ! compute the new mean values

          ! for the stem height (aboveground only)
          values(ipts,ivm,:,iheight) = &
               wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
               ivm,pipe_tune2(ipts,ivm))
         
          ! for the stem diameter
          values(ipts,ivm,:,idiameter) = &
               wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
               ivm,pipe_tune2(ipts,ivm))
         
          ! vertical and horizontal crown diameter
          CALL wood_to_cn_dia(circ_class_biomass(ipts,ivm,:,:,icarbon), &
               circ_class_n(ipts,ivm,:),ivm,values(ipts,ivm,:,icndiahor), &
               values(ipts,ivm,:,icndiaver),pipe_tune2(ipts,ivm))
         
          ! for the crown area
          values(ipts,ivm,:,icnarea) = &
               pi/4*values(ipts,ivm,:,icndiahor)**2
         
          ! for the crown volume
          values(ipts,ivm,:,icnvol) = &
               wood_to_cv(circ_class_biomass(ipts,ivm,:,:,icarbon),&
               circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

       ENDDO ! loop over PFT

    ENDDO ! loop over points

    ! Output tree and crown properties
    ! CCDIAMETER and CCHEIGHT are written in stomate_growth_fun_all
    ! at the same time as some tree ring width variables. They could
    ! as well be written here.
    CALL xios_orchidee_send_field("CCDIAHOR",values(:,:,:,icndiahor))
    CALL xios_orchidee_send_field("CCDIAVER",values(:,:,:,icndiaver))
    CALL xios_orchidee_send_field("CCCROWN",values(:,:,:,icnarea))
    CALL xios_orchidee_send_field("CCCROWNVOL",values(:,:,:,icnvol))

    IF (printlev.GE.3) WRITE(numout,*) 'Leaving derive_biomass_quantities'

 END SUBROUTINE derive_biomass_quantities


END MODULE stomate_stand_structure



