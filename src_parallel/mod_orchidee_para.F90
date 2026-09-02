! ==============================================================================================================================
! MODULE   : mod_orchidee_para
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF      Initialization of MPI and OpenMP parallelization.
!!
!!\n DESCRIPTION  :  This module contains subroutines to be called for the initialization of MPI and OpenMP parallelization. 
!!                   Note that some subroutines are called only for the offline case such as init_orchidee_para and 
!!                   init_orchidee_data_para_driver.
!!
!! SVN              :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parallel/mod_orchidee_para.F90 $
!! $Date: 2025-11-17 09:49:48 +0100 (lun. 17 nov. 2025) $
!! $Revision: 9167 $
!! \n
!_ ================================================================================================================================
MODULE mod_orchidee_para

  USE mod_orchidee_para_var
  USE mod_orchidee_mpi_data
  USE mod_orchidee_omp_data
  USE mod_orchidee_transfert_para
  USE constantes_var, ONLY : printlev
  USE ioipsl_para, ONLY : getin_p

CONTAINS
    
  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_orchidee_para
  !!
  !>\BRIEF	 Initialization of MPI and OpenMP parallelization in offline case
  !!
  !! DESCRIPTION: First subroutine for initialization to be called for the initialization of the MPI and OpenMP parallelization
  !!              in offline mode. This routine will call the successively the initialization for OMP then for MPI.
  !!              We define in this routine the variable "is_root_prc = is_mpi_root AND is_omp_root".
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Init_orchidee_para(communicator, usexios)
    IMPLICIT NONE
    INTEGER,OPTIONAL,INTENT(in) :: communicator	
    LOGICAL,OPTIONAL,INTENT(in) :: usexios
    
    CALL Init_orchidee_omp

    IF ( PRESENT(communicator) ) THEN
       IF ( PRESENT(usexios) ) THEN
          CALL Init_orchidee_mpi(communicator, usexios)
       ELSE
          CALL Init_orchidee_mpi(communicator)
       ENDIF
    ELSE
       CALL Init_orchidee_mpi
    ENDIF

    IF (is_mpi_root .AND. is_omp_root) THEN
       is_root_prc=.TRUE.
    ELSE
       is_root_prc=.FALSE.
    ENDIF
  END SUBROUTINE Init_orchidee_para
    
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_orchidee_data_para_driver
  !!
  !>\BRIEF	 Initialization of variables related to the local domain decomposition called by the offline driver.
  !!
  !! DESCRIPTION: Initialization of variables related to the local domain decomposition.
  !!              This subroutine is only called in offline mode by the driver. 
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Init_orchidee_data_para_driver(nbp,kindex_glo)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbp
    INTEGER,INTENT(IN) :: kindex_glo(nbp)
      
    INTEGER :: first_point
    INTEGER :: last_point
    INTEGER :: nbp_loc
    INTEGER :: nbp_loc_para(0:mpi_size-1)
    INTEGER,ALLOCATABLE :: kindex_loc(:)
    INTEGER :: offset
    INTEGER :: i
    
      
    last_point=0
    
    CALL read_load_balance(nbp,nbp_loc_para)    
    
    DO i=0,mpi_rank
       nbp_loc=nbp_loc_para(i)
       First_point=last_point+1
       Last_point=last_point+nbp_loc
    ENDDO
    
    ALLOCATE(kindex_loc(nbp_loc))
    DO i=1,nbp_loc
       kindex_loc(i)=kindex_glo(i+First_Point-1)
    ENDDO
    
    IF (mpi_rank==0) THEN
       offset=0
    ELSE
       offset=kindex_glo(First_point-1)-MOD(kindex_glo(First_point-1),iim_g)
    ENDIF

    kindex_loc(:)=kindex_loc(:)-offset

    CALL Init_orchidee_data_para(nbp_loc,kindex_loc,offset,omp_size,omp_rank,MPI_COMM_ORCH)
    CALL ipslnlf(new_number=numout)
        
  END SUBROUTINE Init_orchidee_data_para_driver
    
  
  !!  =============================================================================================================================
  !! SUBROUTINE:  Init_orchidee_data_para
  !!
  !>\BRIEF	 Initialization of MPI and OpenMP parallelization.
  !!
  !! DESCRIPTION: Initialization of MPI and OpenMP parallelization.
  !!              This subroutine is called from both the offline driver and from the initialization routine for the coupled mode. 
  !!              This routine will call the successively the initialization for omp and then for mpi.
  !!              We define in this routine the variable "is_root_prc = is_mpi_root AND is_omp_root".
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Init_orchidee_data_para(nbp,kindex,arg_offset,arg_omp_size,arg_omp_rank,COMM)

    IMPLICIT NONE
    INTEGER,INTENT(IN)     :: nbp
    INTEGER,INTENT(IN)     :: kindex(nbp)
    INTEGER,INTENT(IN)     :: arg_offset
    INTEGER,INTENT(IN)     :: arg_omp_size
    INTEGER,INTENT(IN)     :: arg_omp_rank
    INTEGER,INTENT(IN)     :: COMM
    
    INTEGER,SAVE              :: arg_nbp_mpi
    INTEGER,ALLOCATABLE,SAVE  :: kindex_mpi(:)
    LOGICAL                   :: last
    INTEGER                   :: mpi_size
    INTEGER                   :: mpi_rank
    INTEGER                   :: ierr
        
#ifdef CPP_PARA
    CALL MPI_COMM_SIZE(COMM,mpi_size,ierr)
    CALL MPI_COMM_RANK(COMM,mpi_rank,ierr)
#else
    mpi_rank=0
    mpi_size=1
#endif
    
    offset=arg_offset 
    last=.FALSE.
    IF (mpi_rank==mpi_size-1 .AND. arg_omp_rank==arg_omp_size-1) last=.TRUE.
    CALL init_orchidee_omp_data(arg_omp_size,arg_omp_rank,nbp,kindex, offset,last)
    
    IF (is_omp_root) THEN
       arg_nbp_mpi=SUM(nbp_omp_para_nb(:))
       ALLOCATE(kindex_mpi(arg_nbp_mpi))
    ENDIF

    CALL barrier2_omp()
    kindex_mpi(nbp_omp_begin:nbp_omp_end)=kindex(:)+offset
    CALL barrier2_omp()
      
    IF (is_omp_root) THEN      
       kindex_mpi(:)=kindex_mpi(:)-offset
       CALL init_orchidee_mpi_data(arg_nbp_mpi,kindex_mpi,offset,COMM)
       nbp_glo=SUM(nbp_mpi_para(:))
    ENDIF
    CALL barrier2_omp()

    CALL init_gridcells_info

    nbp_loc=nbp

    ! Define is_root_prc
    ! Note that this is already done in init_orchidee_para for the offline case but it is done here again for the coupled case. 
    IF (is_mpi_root .AND. is_omp_root) THEN
       is_root_prc=.TRUE.
    ELSE
       is_root_prc=.FALSE.
    ENDIF
    
    CALL Test_orchidee_para

  END SUBROUTINE Init_orchidee_data_para
    
  !!  =============================================================================================================================
  !! SUBROUTINE:  Set_stdout_file
  !!
  !>\BRIEF	 for each output file will give a unit number for the write function
  !!
  !! DESCRIPTION:	for each output file will give a unit number for the write function
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Set_stdout_file(filename)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: filename
    CHARACTER(len=255) :: fileout
    CHARACTER(len=4)  :: num_mpi
    CHARACTER(len=4)  :: num_omp
    INTEGER,PARAMETER :: base_numout=100
    INTEGER           :: ierr
    LOGICAL           :: printstandout

    !Config Key   = PRINTLEV
    !Config Desc  = Print level for text output
    !Config If    = 
    !Config Help  = Possible values are:
    !Config         0    Nearly no output 
    !Config         1    Minimum writing for long simulations, only master MPI
    !                    and master thread creates an output file. The other
    !                    processors writes to standard output and their level is changed to 0  
    !Config         2    More basic information for long simulations,
    !                    all processors opens an own output text file 
    !Config         3    First debug level, 
    !Config         4    Higher debug level
    !Config Def   = 2
    !Config Units = [0, 1, 2, 3, 4]
    CALL getin_p('PRINTLEV',printlev)

    
    IF (printlev .EQ. 1) THEN
       IF (is_mpi_root.AND.is_omp_root) THEN
          ! Keep output level 1 only for the master processor
          printlev=1
       ELSE
          ! Change output level for all other processors
          printlev=0
       END IF
    END IF

    !Config Key   = PRINTSTANDOUT
    !Config Desc  = Print to standard output
    !Config If    = 
    !Config Help  = Write to standard output instead of opening a text output file out_orchidee 
    !Config Def   = FALSE
    !Config Units = [FLAG]
    printstandout=.FALSE.
    CALL getin_p('PRINTSTANDOUT',printstandout)
    
    IF (is_ok_mpi) THEN
       WRITE(num_mpi,'(I4.4)') mpi_rank
    ENDIF
    
    IF (is_ok_omp) THEN
       WRITE(num_omp,'(I4.4)') omp_rank
    ENDIF
    
     
    IF (is_ok_mpi .AND. is_ok_omp) THEN
       fileout=TRIM(filename)//'_'//num_mpi//'.'//num_omp
       numout=base_numout+omp_rank
    ELSE IF (is_ok_mpi .AND. (.NOT. is_ok_omp)) THEN
       fileout=TRIM(filename)//'_'//num_mpi
       numout=base_numout
    ELSE IF ((.NOT. is_ok_mpi) .AND. is_ok_omp) THEN
       fileout=TRIM(filename)//'_'//num_omp
       numout=base_numout+omp_rank
    ELSE
       fileout=TRIM(filename)
       numout=base_numout
    ENDIF


    IF (printstandout .OR. printlev.EQ.0) THEN
       ! Don't open text output file. All write(numout,*) are sent to the standard output.
       numout=6
    ELSE 
       ! Open text output file. All write(numout,*) are sent to this file.
       OPEN(UNIT=numout,FILE=TRIM(fileout),ACTION='write',STATUS='unknown',FORM='formatted',IOSTAT=ierr) 
       IF (ierr /= 0) THEN
#ifdef CPP_PARA
          CALL MPI_FINALIZE(ierr)
#endif
          WRITE(*,*) "In Set_stdout_file : Erreur can't open file ", filename
          STOP 1
       ENDIF
    END IF
 
    CALL Init_numout_omp(numout)

    IF (printlev>=1) WRITE(numout,*) 'ORCHIDEE start'
  END SUBROUTINE Set_stdout_file
      
      
  !!  =============================================================================================================================
  !! SUBROUTINE:  Test_orchidee_para
  !!
  !>\BRIEF	 
  !!
  !! DESCRIPTION:	
  !!
  !! \n
  !_ ============================================================================================================================== 
  SUBROUTINE Test_orchidee_para

    IMPLICIT NONE

    INTEGER,PARAMETER :: dimsize=3
    REAL :: Array(nbp_loc,dimsize)
    REAL :: Array_glo(nbp_glo,dimsize)
    REAL :: Array_glo_tmp(nbp_glo,dimsize)
    REAL :: Array2D_loc(iim_g,jj_nb)
    REAL :: Array2D_glo(iim_g,jjm_g)
    REAL :: sum1,sum2,sum3
    
    INTEGER :: i,j
    
    DO j=1,dimsize
       DO i=1,nbp_loc
          Array(i,j)=10*j+omp_rank+i*1000
       ENDDO
    ENDDO
      
    CALL gather(Array,Array_glo)
    CALL bcast(Array_glo)
    CALL scatter(Array_glo,array)
    CALL gather(array,array_glo_tmp)
    CALL bcast(array_glo_tmp)    
!    WRITE(*,*) "1) Test parallelism (rank ",mpi_rank,omp_rank,"), Sould be 0 :",SUM(array_glo-array_glo_tmp)

    sum1=SUM(array)
    CALL reduce_sum(sum1,sum2)
    CALL bcast(sum2)
    sum3=SUM(array_glo)
!    WRITE(*,*) "2) Test parallelism (rank ",mpi_rank,omp_rank,"), Sould be 0 :",sum3-sum2
    
    IF (is_omp_root) THEN
       DO j=1,jjm_g
          DO i=1,iim_g
             Array2D_glo(i,j)=(j-1)*iim_g+i
          ENDDO
       ENDDO
        
       array2D_loc(:,:)=0
       CALL scatter2D_mpi(array2D_glo,array2D_loc)
       array2D_glo(:,:)=0
       CALL gather2D_mpi(array2D_loc,array2D_glo)
       CALL bcast_mpi(array2D_glo)
       sum1=SUM(array2D_glo)
       sum2=SUM(array2D_loc)
       CALL reduce_sum_mpi(sum2,sum3)
       CALL bcast_mpi(sum3)
       
!       WRITE(*,*) "3) Test parallelism (rank ",mpi_rank,omp_rank,"), Sould be 0 :",sum3-sum1
    ENDIF
    CALL barrier2_omp()

  END SUBROUTINE  Test_orchidee_para
  
END MODULE mod_orchidee_para
