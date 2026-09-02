! =================================================================================================================================
! MODULE       : structures
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Contains structures and types used throughout the code
!!
!!\n DESCRIPTION: In this module, structural types should be declared.  There are
!!                also routines for allocating and deallocating the memory in structures,
!!                if necessary.
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/structures.f90 $
!! $Date: 2023-07-04 15:56:39 +0200 (mar. 04 juil. 2023) $
!! $Revision: 8071 $
!! \n
!_ ================================================================================================================================

MODULE structures

  USE defprec
  USE constantes
  USE pft_parameters
!-
  IMPLICIT NONE
!-

  ! The purpose of this structure is to contain all
  ! the information needed to fit the effective LAI
  ! as a function of the solar angle.  Doing this
  ! as a structure means that we can change which
  ! function we fit to more easily.  We will have
  ! to allocate this for every grid cell, PFT, and
  ! level.
  INTEGER,PARAMETER :: nparams_laieff = 5 ! The number of parameters in the structure below
  TYPE laieff_type
     REAL(r_std) :: a,b,c,d,e          !! We do a polynomial fit, a+b*x+c*x**2+d*x**3+e*x**4
  END type laieff_type


 CONTAINS

!! ==============================================================================================================================
!! SUBROUTINE   : laieff_type_init
!!
!>\BRIEF        Set the values for all the parameters of this type to zero.
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::laieff_type
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!!
!!
!!
!_ ================================================================================================================================

   SUBROUTINE laieff_type_init(npts, nlevels_loc, laieff_fit)

     !! 0. Variable and parameter declaration

     !! 0.1 Input variables
     INTEGER(i_std), INTENT(in)                                :: npts        !! Domain size - Number of land pixels  (unitless)       
     INTEGER(i_std), INTENT(in)                                :: nlevels_loc     !! The number of vertical levels in the canopy (unitless)    
     !! 0.2 Output variables
     
     !! 0.3 Modified variables
     TYPE(laieff_type),DIMENSION (:,:,:),INTENT(inout)    :: laieff_fit      !! Fitted parameters for the effective LAI
   
     !! 0.4 Local variables

     INTEGER :: ipts,ivm,ilevel

!_ ================================================================================================================================

     DO ipts=1,npts
        DO ivm=1,nvm
           DO ilevel=1,nlevels_loc
              laieff_fit(ipts,ivm,ilevel)%a=zero
              laieff_fit(ipts,ivm,ilevel)%b=zero
              laieff_fit(ipts,ivm,ilevel)%c=zero
              laieff_fit(ipts,ivm,ilevel)%d=zero
              laieff_fit(ipts,ivm,ilevel)%e=zero
           ENDDO
        ENDDO
     ENDDO


   END SUBROUTINE laieff_type_init



 END MODULE structures
