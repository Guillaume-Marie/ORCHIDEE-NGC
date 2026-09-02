! =================================================================================================================================
! MODULE       : dynamic_parameters_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module contains the dynamic parameters.
!!
!!\n DESCRIPTION: 
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_parameters/pft_parameters_var.f90 $
!! $Date: 2023-07-04 15:56:39 +0200 (Tue, 04 Jul 2023) $
!! $Revision: 8071 $
!! \n
!_ ================================================================================================================================

MODULE dynamic_parameters_var

  USE defprec
  
  IMPLICIT NONE

  LOGICAL, SAVE   :: l_first_dynamic_parameters = .TRUE.                 !! To keep first call trace of the module (true/false)
!$OMP THREADPRIVATE(l_first_dynamic_parameters)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: pipe_tune2           !! height=pipe_tune2 * diameter**pipe_tune3
!$OMP THREADPRIVATE(pipe_tune2)      

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: alpha_self_thinning  !! Coefficient of the self-thinning relationship D=alpha*N^beta 
                                                                         !! estimated from German, French, Spanish and Swedish 
                                                                         !! forest inventories
!$OMP THREADPRIVATE(alpha_self_thinning)

END MODULE dynamic_parameters_var
