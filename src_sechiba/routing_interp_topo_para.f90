MODULE routing_interp_topo_para
#ifdef XIOS
 PRIVATE

 INTEGER,SAVE  :: ni
 !$OMP THREADPRIVATE(ni)
 INTEGER,SAVE  :: nj
 !$OMP THREADPRIVATE(nj)

 INTEGER,SAVE  :: ni_glo
 !$OMP THREADPRIVATE(ni_glo)
 INTEGER,SAVE  :: nj_glo
 !$OMP THREADPRIVATE(nj_glo)
 INTEGER,SAVE  :: nbpt_r 
 !$OMP THREADPRIVATE(nbpt_r)
 INTEGER,SAVE  :: nbpt_rp1 
 !$OMP THREADPRIVATE(nbpt_rp1)
 INTEGER,SAVE,ALLOCATABLE :: r_to_rp1(:)
 !$OMP THREADPRIVATE(r_to_rp1)
 INTEGER,SAVE,ALLOCATABLE :: rp1_to_r(:)
 !$OMP THREADPRIVATE(rp1_to_r)

 TYPE mpi_neighbours
   INTEGER             :: rank
   INTEGER             :: n
   INTEGER,ALLOCATABLE :: index(:)
   INTEGER,ALLOCATABLE :: ibuffer(:)
   REAL,ALLOCATABLE :: rbuffer(:)
 END TYPE
 
 TYPE(mpi_neighbours),SAVE, ALLOCATABLE :: send_neighbours(:)
 TYPE(mpi_neighbours),SAVE, ALLOCATABLE :: recv_neighbours(:)
 
 INTERFACE update_halo
   PROCEDURE update_halo_integer, update_halo_real
 END INTERFACE

 PUBLIC routing_para_initialize, update_halo, r_to_rp1, rp1_to_r

CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_para_initialize
  !!
  !>\BRIEF	         Initialize the routing_interp_topo_para module
  !!
  !! DESCRIPTION:        Initialize the routing_interp_topo_para module.
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE routing_para_initialize
  USE mod_orchidee_para, ONLY : mpi_rank,mpi_size, MPI_COMM_ORCH, MPI_INT_ORCH, is_omp_root
  USE xios
  IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL, ALLOCATABLE :: my_rank(:)
    REAL, ALLOCATABLE :: global_ranks(:)
    REAL, ALLOCATABLE :: my_index(:)
    REAL, ALLOCATABLE :: global_index(:)
    LOGICAL, ALLOCATABLE :: inner(:)      ! is the pixel (from the domain local+halo) in the local domain (not halo)?
    INTEGER :: nb_send_neighbours
    INTEGER :: nb_recv_neighbours
    INTEGER :: nb_points_to_recv
    INTEGER, ALLOCATABLE :: list_rank_to_recv(:)
     INTEGER :: i,j,k,ij,ind,r,rp1
    INTEGER :: rank, pos_rank
    INTEGER,ALLOCATABLE :: send_request(:)
    INTEGER,ALLOCATABLE :: recv_request(:)
    INTEGER,ALLOCATABLE :: send_status(:,:)
    INTEGER,ALLOCATABLE :: recv_status(:,:)
    INTEGER :: nb_points_ranks(0:mpi_size-1)
    INTEGER :: ierr
    INTEGER :: test

    IF (is_omp_root) THEN
      CALL xios_get_domain_attr("routing_domain", ni=ni, nj=nj, ni_glo=ni_glo, nj_glo=nj_glo)    ! get routing domain dimension
  
      nbpt_r= ni*nj                                                
      nbpt_rp1= (ni+2)*(nj+2)
      
      ALLOCATE(my_rank(nbpt_r))
      ALLOCATE(global_ranks(nbpt_rp1))
      ALLOCATE(my_index(nbpt_r))
      ALLOCATE(global_index(nbpt_rp1))
      ALLOCATE(r_to_rp1(nbpt_r))
      ALLOCATE(rp1_to_r(nbpt_rp1))
      ALLOCATE(inner(nbpt_rp1))
  
      my_rank(:)=mpi_rank
      rp1_to_r(:)=-1
      inner(:)=.FALSE.
  
      DO j=1,nj
        DO i=1,ni
          r=((j-1)*ni)+i
          rp1=(j+1-1)*(ni+2)+i+1 
          r_to_rp1(r)=rp1
          rp1_to_r(rp1)=r
          my_index(r) = rp1
          inner(rp1)=.TRUE.
        ENDDO
      ENDDO
  
      CALL xios_send_field("my_rank",my_rank)
      CALL xios_recv_field("global_ranks",global_ranks)
      CALL xios_send_field("my_index",my_index)
      CALL xios_recv_field("global_index",global_index)
  
      ! number of point received by rank 
      nb_points_ranks(:)=0
      DO rp1=1,nbpt_rp1
        nb_points_ranks(global_ranks(rp1)) = nb_points_ranks(global_ranks(rp1))+1
      ENDDO
      
      nb_recv_neighbours=0
      nb_points_to_recv=0
      DO rank=0, mpi_size-1
        IF (nb_points_ranks(rank)/=0) THEN
          nb_recv_neighbours=nb_recv_neighbours+1
          nb_points_to_recv = nb_points_to_recv + nb_points_ranks(rank)
        ENDIF
      ENDDO
      
      ALLOCATE(recv_neighbours(nb_recv_neighbours))
      ALLOCATE(list_rank_to_recv(0:mpi_size-1))
      k=1
      DO rank=0, mpi_size-1
        IF (nb_points_ranks(rank)/=0) THEN
          recv_neighbours(k)%rank = rank
          recv_neighbours(k)%n = 0
          ALLOCATE(recv_neighbours(k)%index(nb_points_ranks(rank)))
          ALLOCATE(recv_neighbours(k)%ibuffer(nb_points_ranks(rank)))
          ALLOCATE(recv_neighbours(k)%rbuffer(nb_points_ranks(rank)))
          list_rank_to_recv(rank) = k
          k=k+1
        ENDIF
      ENDDO    
      
      DO rp1=1,nbpt_rp1
        rank = global_ranks(rp1)
        IF (.NOT. inner(rp1)) THEN
          k = list_rank_to_recv(rank) 
          recv_neighbours(k)%n=recv_neighbours(k)%n+1
          recv_neighbours(k)%index(recv_neighbours(k)%n) = rp1
        ENDIF
      ENDDO
  
  
      nb_points_ranks(:) = 0
      DO k=1,size(recv_neighbours)
        nb_points_ranks(recv_neighbours(k)%rank) = 1
      ENDDO
      
      CALL MPI_ALLREDUCE(MPI_IN_PLACE , nb_points_ranks, mpi_size, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
  
      nb_send_neighbours = nb_points_ranks(mpi_rank)
      ALLOCATE(send_neighbours(nb_send_neighbours))
  
  
      ALLOCATE(send_request(size(recv_neighbours)))
      ALLOCATE(send_status(MPI_STATUS_SIZE,size(recv_neighbours)))
      ALLOCATE(recv_request(size(send_neighbours)))
      ALLOCATE(recv_status(MPI_STATUS_SIZE, size(send_neighbours)))
      

      DO k=1,size(recv_neighbours)
        CALL MPI_ISEND(recv_neighbours(k)%n, 1, MPI_INT_ORCH, recv_neighbours(k)%rank, 0, MPI_COMM_ORCH, send_request(k) , ierr )
      ENDDO

      DO k=1,size(send_neighbours)
        CALL MPI_IRecv(send_neighbours(k)%n, 1, MPI_INT_ORCH, MPI_ANY_SOURCE, 0, MPI_COMM_ORCH, recv_request(k) , ierr )
      ENDDO

      CALL MPI_WAITALL(size(recv_neighbours),send_request, send_status, ierr)
      CALL MPI_WAITALL(size(send_neighbours),recv_request, recv_status, ierr)
    
      DO k=1,size(send_neighbours)
        send_neighbours(k)%rank = recv_status(MPI_SOURCE, k)
        ALLOCATE(send_neighbours(k)%index(send_neighbours(k)%n))
        ALLOCATE(send_neighbours(k)%ibuffer(send_neighbours(k)%n))
        ALLOCATE(send_neighbours(k)%rbuffer(send_neighbours(k)%n))
      ENDDO

      DO k=1,size(recv_neighbours)
        DO ij=1, recv_neighbours(k)%n
          ind = recv_neighbours(k)%index(ij)
          recv_neighbours(k)%ibuffer(ij) = global_index(ind)
        ENDDO
      ENDDO

      DO k=1,size(recv_neighbours)
        CALL MPI_ISEND(recv_neighbours(k)%ibuffer, recv_neighbours(k)%n, MPI_INT_ORCH, recv_neighbours(k)%rank, 1, MPI_COMM_ORCH, send_request(k) , ierr )
      ENDDO

      DO k=1,size(send_neighbours)
        CALL MPI_IRECV(send_neighbours(k)%index, send_neighbours(k)%n, MPI_INT_ORCH, send_neighbours(k)%rank, 1, MPI_COMM_ORCH, recv_request(k) , ierr )
      ENDDO

      CALL MPI_WAITALL(size(recv_neighbours),send_request, send_status, ierr)
      CALL MPI_WAITALL(size(send_neighbours),recv_request, recv_status, ierr)

      global_ranks(:)=1000000000
      DO r=1,nbpt_r
        global_ranks(r_to_rp1(r)) = mpi_rank 
      ENDDO
    
      CALL update_halo_real(global_ranks)
      
    ENDIF ! is_omp_root

  END SUBROUTINE routing_para_initialize


  !!  =============================================================================================================================
  !! SUBROUTINE:         update_halo_integer
  !!
  !>\BRIEF	         Procedure for integers of interface update_halo
  !!
  !! DESCRIPTION:        Procedure for integers of interface update_halo. Must be called from an omp_root thread.
  !!                     
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE update_halo_integer(field)
  USE mod_orchidee_para, ONLY :  MPI_COMM_ORCH, MPI_INT_ORCH
  IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(INOUT) :: field(nbpt_rp1)
    INTEGER :: send_request(size(send_neighbours))
    INTEGER :: send_status(MPI_STATUS_SIZE,size(send_neighbours))
    INTEGER :: recv_request(size(recv_neighbours))
    INTEGER :: recv_status(MPI_STATUS_SIZE, size(recv_neighbours))
    INTEGER :: ij, k,ind
    INTEGER :: ierr
   
    DO k=1,size(send_neighbours)
      DO ij=1, send_neighbours(k)%n
        ind = send_neighbours(k)%index(ij)
        send_neighbours(k)%ibuffer(ij) = field(ind)
      ENDDO
    ENDDO

    DO ij=1,size(send_neighbours)
      CALL MPI_ISEND(send_neighbours(ij)%ibuffer, send_neighbours(ij)%n, MPI_INT_ORCH, send_neighbours(ij)%rank, 0, MPI_COMM_ORCH, send_request(ij) , ierr )
    ENDDO

    DO ij=1,size(recv_neighbours)
      CALL MPI_IRECV(recv_neighbours(ij)%ibuffer, recv_neighbours(ij)%n, MPI_INT_ORCH, recv_neighbours(ij)%rank, 0, MPI_COMM_ORCH, recv_request(ij) , ierr )
    ENDDO

    CALL MPI_WAITALL(size(send_neighbours),send_request, send_status, ierr)
    CALL MPI_WAITALL(size(recv_neighbours),recv_request, recv_status, ierr)

    DO k=1,size(recv_neighbours)
      DO ij=1, recv_neighbours(k)%n
        ind = recv_neighbours(k)%index(ij)
        field(ind) = recv_neighbours(k)%ibuffer(ij) 
      ENDDO
    ENDDO

  END SUBROUTINE update_halo_integer

  !!  =============================================================================================================================
  !! SUBROUTINE:         update_halo_real
  !!
  !>\BRIEF	         Procedure for reals of interface update_halo
  !!
  !! DESCRIPTION:        Procedure for reals of interface update_halo. Must be called from an omp_root thread.
  !!                     
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE update_halo_real(field)
  USE mod_orchidee_para, ONLY :  MPI_COMM_ORCH, MPI_REAL_ORCH
  IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL, INTENT(INOUT) :: field(nbpt_rp1)
    INTEGER :: send_request(size(send_neighbours))
    INTEGER :: send_status(MPI_STATUS_SIZE,size(send_neighbours))
    INTEGER :: recv_request(size(recv_neighbours))
    INTEGER :: recv_status(MPI_STATUS_SIZE, size(recv_neighbours))
    INTEGER :: ij, k,ind
    INTEGER :: ierr

    DO k=1,size(send_neighbours)
      DO ij=1, send_neighbours(k)%n
        ind = send_neighbours(k)%index(ij)
        send_neighbours(k)%rbuffer(ij) = field(ind)
      ENDDO
    ENDDO
   
    DO ij=1,size(send_neighbours)
      CALL MPI_ISEND(send_neighbours(ij)%rbuffer, send_neighbours(ij)%n, MPI_REAL_ORCH, send_neighbours(ij)%rank, 0, MPI_COMM_ORCH, send_request(ij) , ierr )
    ENDDO

    DO ij=1,size(recv_neighbours)
      CALL MPI_IRECV(recv_neighbours(ij)%rbuffer, recv_neighbours(ij)%n, MPI_REAL_ORCH, recv_neighbours(ij)%rank, 0, MPI_COMM_ORCH, recv_request(ij) , ierr )
    ENDDO

    CALL MPI_WAITALL(size(send_neighbours),send_request, send_status, ierr)
    CALL MPI_WAITALL(size(recv_neighbours),recv_request, recv_status, ierr)

    DO k=1,size(recv_neighbours)
      DO ij=1, recv_neighbours(k)%n
        ind = recv_neighbours(k)%index(ij)
        field(ind) = recv_neighbours(k)%rbuffer(ij) 
      ENDDO
    ENDDO

  END SUBROUTINE update_halo_real
#endif
END MODULE routing_interp_topo_para
