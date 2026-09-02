! =================================================================================================================================
! MODULE       : utils
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Modules containing divers gerneric functions and subroutines
!!
!!\n DESCRIPTION: Modules containing divers gerneric functions and subroutines. Contains subroutine:
!!                - nccheck for handeling netcdf output message.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_global/utils.f90 $
!! $Date: 2023-07-04 15:51:03 +0200 (mar. 04 juil. 2023) $
!! $Revision: 8070 $
!_
!================================================================================================================================


MODULE utils 

  USE netcdf
  USE defprec
  USE ioipsl_para

  IMPLICIT NONE

  PRIVATE
  PUBLIC nccheck 

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: nccheck 
!!
!>\BRIEF        Check for netcdf exit status 
!!
!! DESCRIPTION  : Launch an orchidee error message if status variable contains a netcdf error
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE nccheck(status)
    INTEGER(i_std), INTENT (IN)         :: status
    CHARACTER(LEN=200)                  :: mesg
    
    IF(status /= nf90_noerr) THEN
      
      WRITE(numout, *) trim(nf90_strerror(status))
      CALL ipslerr_p(3, 'nccheck', 'Netcdf error', 'Check out_orchide_XXXX output files', 'for more information')
    END IF  
  END SUBROUTINE nccheck

END MODULE utils 
