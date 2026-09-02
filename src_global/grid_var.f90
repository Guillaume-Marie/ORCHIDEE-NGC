! ===============================================================================================================================
! MODULE       : grid_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF
!!
!! \n DESCRIPTION : This module define variables for the grid module. 
!!                The module is already USE in module grid. Therefor no need to use it seperatly if use grid is already done.

!!
!! RECENT CHANGE(S): These variables were previously in grid module. They have been moved here to avoid dependency 
!!                   problems when the variables are needed in the parallelization modules. 
!!
!!
!! SVN
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_global/grid_var.f90 $
!! $Date: 2023-07-04 16:36:22 +0200 (mar. 04 juil. 2023) $
!! $Revision: 8073 $
!_ ================================================================================================================================

MODULE grid_var

  USE defprec
  IMPLICIT NONE

  !=================================================================================
  ! Horizontal grid information
  !=================================================================================
  INTEGER,PARAMETER                                  :: unstructured=0   !! Index for unstructured grid
  INTEGER,PARAMETER                                  :: regular_lonlat=1 !! Index for regular longitude-latitude grid
  INTEGER,PARAMETER                                  :: regular_xy=2     !! Index for regular grid projected on X-Y
  INTEGER,SAVE                                       :: grid_type        !! grid type : unstructured, regular_lonlat or regular_xy 
!$OMP THREADPRIVATE(grid_type)
  INTEGER(i_std), SAVE                               :: ncells           !! numbers of cells (including masked cells)
!$OMP THREADPRIVATE(ncells)
  INTEGER(i_std), SAVE                               :: nvertex          !! numbers of vertex by cells
!$OMP THREADPRIVATE(nvertex)
  INTEGER, SAVE                                      :: NbSegments       !! Number of segments in the polygone defining the grid box, same as nvertex
!$OMP THREADPRIVATE(NbSegments)
  INTEGER, SAVE                                      :: NbNeighb         !! Number of neighbours for each grid box
!$OMP THREADPRIVATE(NbNeighb)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE       :: longitude        !! longitude for all cells (including masked cells)
!$OMP THREADPRIVATE(longitude)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE       :: latitude         !! latitude for all cells (including masked cells)
!$OMP THREADPRIVATE(latitude)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE     :: bounds_lon       !! longitude boundairies for all cells (including masked cells)
!$OMP THREADPRIVATE(bounds_lon)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE     :: bounds_lat       !! latitude boundairies for all cells (including masked cells)
!$OMP THREADPRIVATE(bounds_lat)
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:), SAVE    :: ind_cell_glo     !! Cells order give by indice form global grid reference
!$OMP THREADPRIVATE(ind_cell_glo)


END MODULE grid_var
