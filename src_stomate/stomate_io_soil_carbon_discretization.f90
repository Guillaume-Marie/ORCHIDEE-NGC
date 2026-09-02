!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_io_soil_carbon_discretization.f90 $ 
!! $Date: 2026-04-03 18:22:57 +0200 (ven. 03 avril 2026) $
!! $Revision: 9463 $
! IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
MODULE stomate_io_soil_carbon_discretization
  !---------------------------------------------------------------------
  !-
  !- 
  !-
  !---------------------------------------------------------------------
  USE netcdf
  USE defprec
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE vertical_soil
  USE mod_orchidee_para
  USE ioipsl_para 
  USE utils      ! nccheck
#ifdef CPP_PARA
  USE mpi
#endif
  !-
  IMPLICIT NONE
  !-
  PRIVATE
  PUBLIC stomate_io_soil_carbon_discretization_write, stomate_io_soil_carbon_discretization_read
  !-
  ! TO CHECK, stomate_finalize also uses this type of var
  INTEGER,PARAMETER                              :: r_typ = NF90_REAL8   !! Specify data format (server dependent)
  !-
CONTAINS
  !-
  !===
  !-
  !! ================================================================================================================================
  !! SUBROUTINE 	: stomate_io_soil_carbon_discretization_write 
  !!
  !>\BRIEF        Writes stomate permafrost carbon data into a netcdf file
  !!
  !! DESCRIPTION  : It writes into a netcdf files in parallel mode all necessary
  !!                 variables required for spinup (forecesoil)
  !!                
  !!                
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE stomate_io_soil_carbon_discretization_write (Cforcing_permafrost_name, & 
       nbp_glo,            start_px,           length_px,      nparan,     nbyear, &
       index_g,                                      &
       clay,               depth_organic_soil, lalo,                               &
       snowdz_2pfcforcing, snowrho_2pfcforcing, som_input_2pfcforcing, &
       tsurf_2pfcforcing,  pb_2pfcforcing,     snow_2pfcforcing, &
       tprof_2pfcforcing,  fbact_2pfcforcing,  veget_max_2pfcforcing, &
       rprof_2pfcforcing,  hslong_2pfcforcing, CN_target_2pfcforcing, &
       n_mineralisation_2pfcforcing, &
       wtp_pt)


    CHARACTER(LEN=100), INTENT(in)              :: Cforcing_permafrost_name !! Name of permafrost forcing file
    INTEGER(i_std), INTENT(in)                  :: nbp_glo !nbp_glo is the number of global continental points
    INTEGER(i_std), INTENT(in)                  :: start_px ! Start land point/pixex respect to nbp_glo
    INTEGER(i_std), INTENT(in)                  :: length_px ! Length of lands point/pixel to write
    INTEGER(i_std), INTENT(in)                  :: nparan ! Number of forcesoil timesteps  
    INTEGER(i_std), INTENT(in)                  :: nbyear ! Number of years saved for carbon spinup 
    INTEGER(i_std),DIMENSION(:),INTENT(in)      :: index_g             !! Indices of the terrestrial pixels only (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)       :: clay                   !! Clay fraction of soil (0-1, unitless), parallel 
    REAL(r_std), DIMENSION(:), INTENT(in)       :: depth_organic_soil !! Depth at which there is still organic matter (m)!
    REAL(r_std), DIMENSION(:,:),INTENT(in)      :: lalo              !! Geographical coordinates (latitude,longitude) 
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: snowdz_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: snowrho_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(in) :: som_input_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in )     :: tsurf_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in)      :: pb_2pfcforcing
    REAL(r_std),DIMENSION(:,:), INTENT(in)      :: snow_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: tprof_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: fbact_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: hslong_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: veget_max_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: rprof_2pfcforcing
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(in)  :: CN_target_2pfcforcing !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1)   
    REAL(r_std),DIMENSION(:,:,:), INTENT(in)    :: n_mineralisation_2pfcforcing !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
    REAL(r_std),DIMENSION(:,:), INTENT(in)      :: wtp_pt                      !! water table above surface

    ! Local Variables
    INTEGER(i_std)                              :: ier, n_directions, i
    INTEGER(i_std)                              :: start(1), ncount(1), start_2d(2), ncount_2d(2), inival, endval 
    INTEGER(i_std)                              :: start_4d(4), ncount_4d(4), start_3d(3), ncount_3d(3)
    INTEGER(i_std)                              :: start_5d(5), ncount_5d(5)
    INTEGER(i_std),DIMENSION(10)                :: d_id                     !! List each netcdf dimension
    INTEGER(i_std)                              :: vid                      !! Variable identifer of netCDF (unitless)
    INTEGER(i_std)                              :: Cforcing_permafrost_id   !! Permafrost file identifer 

    ! Create file 
#ifdef CPP_PARA
    ier = NF90_CREATE (TRIM(Cforcing_permafrost_name),IOR(NF90_NETCDF4,NF90_MPIIO), &
         Cforcing_permafrost_id, comm=MPI_COMM_ORCH, info=MPI_INFO_NULL)
#else
    ier = NF90_CREATE (TRIM(Cforcing_permafrost_name),NF90_NETCDF4, &
         Cforcing_permafrost_id)
#endif
    IF (ier /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'stomate_finalize', &
            &        'PROBLEM creating Cforcing_permafrost file', &
            &        NF90_STRERROR(ier),'')
    END IF


    ! Add variable attribute
    ! Note ::nbp_glo is the number of global continental points
    CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
         &                           'kjpindex',REAL(nbp_glo,r_std)))
    CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
         &                           'nparan',REAL(nparan,r_std)))
    CALL nccheck( NF90_PUT_ATT (Cforcing_permafrost_id,NF90_GLOBAL, &
         &                           'nbyear',REAL(nbyear,r_std)))

    ! Add new dimension, variables values from USE
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'points',nbp_glo,d_id(1)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'carbtype',ncarb,d_id(2)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'vegtype',nvm,d_id(3)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'level',ngrnd,d_id(4)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'time_step',NF90_UNLIMITED,d_id(5)))
    n_directions=2
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'direction',n_directions,d_id(6)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'elements',nelements,d_id(7)))
    CALL nccheck( NF90_DEF_DIM (Cforcing_permafrost_id,'snowlevel',nsnow,d_id(8)))

    ! Add new variable
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'points', r_typ,d_id(1),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'carbtype', r_typ,d_id(2),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'vegtype', r_typ,d_id(3),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'level', r_typ,d_id(4),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'time_step',r_typ,d_id(5),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'direction',r_typ,d_id(6),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'elements',r_typ,d_id(7),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'snowlevel', r_typ,d_id(8),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'index',r_typ,d_id(1),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'clay',r_typ,d_id(1),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'depth_organic_soil',r_typ,d_id(1),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'lalo',      r_typ, &
         (/ d_id(1), d_id(6) /),vid))
    !--time-invariant
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'zz_deep',r_typ,d_id(4),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'zz_coef_deep',r_typ,d_id(4),vid))
    !
    IF (ok_peat_NoDiscretisation) THEN
       CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'wtp_pt',r_typ, &
            &                        (/ d_id(1),d_id(5)/),vid))
    ENDIF
    !
    !--3layers snow 
    CALL nccheck( NF90_DEF_VAR(Cforcing_permafrost_id,'snowdz',r_typ,(/ d_id(1),d_id(8),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR(Cforcing_permafrost_id,'snowrho',r_typ,(/ d_id(1),d_id(8),d_id(5) /),vid))
    !--time-varying 
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'som_input',r_typ, &
         &                        (/ d_id(1),d_id(2),d_id(3),d_id(7), d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'pb',r_typ, & 
         &                        (/ d_id(1),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'snow',r_typ, &
         &                        (/ d_id(1),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'tprof',r_typ, &
         &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'fbact',r_typ, &
         &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'hslong',r_typ, &
         &                        (/ d_id(1),d_id(4),d_id(3),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'veget_max',r_typ, &
         &                        (/ d_id(1),d_id(3),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'rprof',r_typ, &
         &                        (/ d_id(1),d_id(3),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'tsurf',r_typ, &
         &                        (/ d_id(1),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'CN_target',r_typ, &
         &                        (/ d_id(1),d_id(3), d_id(2),d_id(5) /),vid))
    CALL nccheck( NF90_DEF_VAR (Cforcing_permafrost_id,'n_mineralisation',r_typ, &
         &                        (/ d_id(1),d_id(3), d_id(5) /),vid))
    CALL nccheck( NF90_ENDDEF (Cforcing_permafrost_id))

    ! Write data
    start=(/ start_px /)
    ncount=(/ length_px /)
    inival=start_px
    endval=start_px + length_px !length_px_end(mpi_rank)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'points',vid) )
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &  (/(REAL(i,r_std),i=inival,endval)/), &
         &  start=start, count=ncount) )

    ! no point to make parallel calls
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'carbtype',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &                        (/(REAL(i,r_std),i=1,ncarb)/)))
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'vegtype',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &                            (/(REAL(i,r_std),i=1,nvm)/)))
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'level',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &                        (/(REAL(i,r_std),i=1,ngrnd)/)))
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'time_step',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &                       (/(REAL(i,r_std),i=1,nparan*nbyear)/)))
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'index',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, REAL(index_g,r_std) ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'zz_deep',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, znt ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'zz_coef_deep',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, zlt ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowlevel',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, &
         &                        (/(REAL(i,r_std),i=1,nsnow)/)))
    ! Parallel writes
    start=(/ start_px /)
    ncount=(/ length_px /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'clay',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, clay, start=start, count=ncount  ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'depth_organic_soil',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, depth_organic_soil, start=start, count=ncount  ))

    start_2d=(/ start_px,1 /)
    ncount_2d=(/ length_px,2 /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'lalo',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, lalo, start=start_2d, count=ncount_2d ))

    !peatland++
    start_2d=(/ start_px,1 /)
    ncount_2d=(/ length_px,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID(Cforcing_permafrost_id,'wtp_pt',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, wtp_pt, start=start_2d ,count=ncount_2d ))

    ! putting 3 snow layers
    start_3d=(/ start_px,1,1 /)
    ncount_3d=(/ length_px,nsnow,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID(Cforcing_permafrost_id,'snowdz',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snowdz_2pfcforcing, start=start_3d ,count=ncount_3d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowrho',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snowrho_2pfcforcing, &
         start=start_3d ,count=ncount_3d ))

    start_5d=(/ start_px,1,1,1,1 /)
    ncount_5d=(/ length_px,ncarb,nvm,nelements,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'som_input',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,som_input_2pfcforcing, &
         start=start_5d ,count=ncount_5d ))

    start_2d=(/ start_px,1 /)
    ncount_2d=(/ length_px,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tsurf',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, tsurf_2pfcforcing, start=start_2d, count=ncount_2d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'pb',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, pb_2pfcforcing, start=start_2d, count=ncount_2d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snow',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, snow_2pfcforcing, start=start_2d, count=ncount_2d ))

    start_4d=(/ start_px,1,1,1 /)
    ncount_4d=(/ length_px,ngrnd,nvm,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tprof',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, tprof_2pfcforcing ,start=start_4d, count=ncount_4d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'fbact',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, fbact_2pfcforcing,start=start_4d, count=ncount_4d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'hslong',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, hslong_2pfcforcing,start=start_4d, count=ncount_4d ))

    start_3d=(/ start_px,1,1 /)
    ncount_3d=(/ length_px,nvm,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'veget_max',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,veget_max_2pfcforcing, start=start_3d, count=ncount_3d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'rprof',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, rprof_2pfcforcing, start=start_3d, count=ncount_3d ))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'n_mineralisation',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid, n_mineralisation_2pfcforcing, start=start_3d, count=ncount_3d ))

    start_4d=(/ start_px,1,1,1 /)
    ncount_4d=(/ length_px,nvm,ncarb,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'CN_target',vid))
    CALL nccheck( NF90_VAR_PAR_ACCESS(Cforcing_permafrost_id, vid, NF90_COLLECTIVE))
    CALL nccheck( NF90_PUT_VAR (Cforcing_permafrost_id,vid,CN_target_2pfcforcing, start=start_4d, count=ncount_4d ))


    ! Finish netcdf file management
    CALL nccheck( NF90_CLOSE (Cforcing_permafrost_id) )

  END SUBROUTINE stomate_io_soil_carbon_discretization_write
  !-
  !===
  !-
  SUBROUTINE stomate_io_soil_carbon_discretization_read (Cforcing_permafrost_name,   &
       nparan,             nbyear,     start_px,   length_px,      &
       som_input,   pb,         snow,       tsurf,          & 
       tprof,              fbact,      hslong,     rprof,          &
       lalo,               snowdz,     snowrho,    veget_max, &
       CN_target, n_mineralisation )

    ! Input Variables
    CHARACTER(LEN=100), INTENT(in)               :: Cforcing_permafrost_name !! Name of permafrost forcing file
    INTEGER(i_std), INTENT(in)                   :: start_px ! Start land point/pixex respect to nbp_glo
    INTEGER(i_std), INTENT(in)                   :: length_px ! Length of lands point/pixel to write
    INTEGER(i_std), INTENT(in)                   :: nparan, nbyear 

    ! Output variables
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(out):: som_input
    REAL(r_std), DIMENSION(:,:), INTENT(out)     :: pb 
    REAL(r_std), DIMENSION(:,:), INTENT(out)     :: snow
    REAL(r_std), DIMENSION(:,:), INTENT(out)     :: tsurf
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: tprof
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: fbact 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out) :: hslong
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: rprof
    REAL(r_std), DIMENSION(:,:), INTENT(out)     :: lalo
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: snowdz
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: snowrho 
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)   :: veget_max 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)   :: CN_target !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1)   
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)  :: n_mineralisation

    ! Local Variables
    INTEGER(i_std)                               :: start_2d(2), count_2d(2) 
    INTEGER(i_std)                               :: start_4d(4), count_4d(4), start_3d(3), count_3d(3)
    INTEGER(i_std)                               :: start_5d(5), count_5d(5)
    INTEGER(i_std)                               :: v_id                      !! Variable identifer of netCDF (unitless)
    INTEGER(i_std)                               :: Cforcing_permafrost_id   !! Permafrost file identifer 
    !-
    ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
    ! and allocate variables.
    !-
#ifdef CPP_PARA
    CALL nccheck( NF90_OPEN (TRIM(Cforcing_permafrost_name),IOR(NF90_NOWRITE, NF90_MPIIO),Cforcing_permafrost_id, &
         & comm = MPI_COMM_ORCH, info = MPI_INFO_NULL ))
#else
    CALL nccheck( NF90_OPEN (TRIM(Cforcing_permafrost_name),NF90_NOWRITE,Cforcing_permafrost_id))
#endif

    start_5d = (/ start_px, 1, 1, 1, 1 /)
    count_5d = (/ length_px, ncarb, nvm, nelements, nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'som_input',v_id))
    CALL nccheck( NF90_GET_VAR   (Cforcing_permafrost_id,v_id,som_input,  &
         &  start = start_5d, count = count_5d ))

    start_2d=(/ start_px, 1 /)
    count_2d=(/ length_px, nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'pb',v_id ))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,pb, &
         & start=start_2d, count=count_2d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snow',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snow, &
         & start=start_2d, count=count_2d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tsurf',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,tsurf, &
         & start=start_2d, count=count_2d))

    start_4d=(/ start_px,1,1,1 /)
    count_4d=(/ length_px,ngrnd,nvm,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'tprof',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,tprof, &
         & start=start_4d, count=count_4d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'fbact',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,fbact, &
         & start=start_4d, count=count_4d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'hslong',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,hslong, &
         & start=start_4d, count=count_4d))

    start_3d=(/ start_px,1,1 /)
    count_3d=(/ length_px,nvm,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'veget_max',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,veget_max, &
         & start=start_3d, count=count_3d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'rprof',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,rprof, &
         & start=start_3d, count=count_3d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'n_mineralisation',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id, n_mineralisation, &
         & start=start_3d, count=count_3d))

    start_2d=(/ start_px, 1 /)
    count_2d=(/ length_px, 2 /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'lalo',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,lalo, &
         & start=start_2d, count=count_2d))

    start_3d=(/ start_px,1,1 /)
    count_3d=(/ length_px,nsnow,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowdz',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snowdz, &
         & start=start_3d, count=count_3d))

    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'snowrho',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,snowrho, &
         & start=start_3d, count=count_3d))

    start_4d=(/ start_px,1,1,1 /)
    count_4d=(/ length_px,nvm,ncarb,nparan*nbyear /)
    CALL nccheck( NF90_INQ_VARID (Cforcing_permafrost_id,'CN_target',v_id))
    CALL nccheck( NF90_GET_VAR (Cforcing_permafrost_id,v_id,CN_target, &
         & start=start_4d, count=count_4d))
    !- Close Netcdf carbon permafrost file reference
    CALL nccheck( NF90_CLOSE (Cforcing_permafrost_id))

  END SUBROUTINE stomate_io_soil_carbon_discretization_read

END MODULE stomate_io_soil_carbon_discretization
