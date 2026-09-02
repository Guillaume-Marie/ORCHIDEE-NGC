!  ==============================================================================================================================\n
!  MODULE forcing_tools	: This module concentrates on the temporal interpolation of the forcing for ORCHIDEE. 
!                         It provides basic service for the grid when this is provided in the forcing file. The main
!                         work for the grid is done in glogrid.f90. The approach of forcing_tools to handle the time
!                         aspect of the forcing is to read as many time steps as possible in memory and then
!                         interpolate that to the time interval requested by the calling program.
!                         The data is read on root_proc but then distributed over all processors according to the
!                         domain decomposition of ORCHIDEE. This allows to be more efficient in the temporal interpolation.
!                         It is important to keep in mind that forcing_tools works on time intervals. So the request for data 
!                         of ORCHIDEE as to be for an interval and the forcing file needs to have a description of the time interval
!                         over which the forcing is valid.
!                         The general description of how the attributes needed in the netCDF file for describing the cell_methods
!                         for time are provided in this document :
!                          https://forge.ipsl.jussieu.fr/orchidee/attachment/wiki/Documentation/Forcings/Description_Forcing_Files.pdf
!
!                         The most important routines of foring_tools are forcing_open and forcing_getvalues
!
!                       forcing_integration_time : Computes the interval over which the simulation should be carried out.
!                       forcing_open : Opens the forcing files and extract the main information.
!                       forcing_getvalues : Gets the forcing data for a time interval.
!                       forcing_close : Closes the forcing file
!                       forcing_printdate : A tool to print the dates in human readable form.
!                       forcing_printpoint : Print the values for a given point in time.
!                       forcing_givegridsize : Allows other programs to get the dimensions of the forcing grid.
!                       forcing_getglogrid : Allows other programs to get the spatial grid of the forcing.
!                       forcing_givegrid : Returns the description of the grid.
!                       forcing_zoomgrid : Extract a sub-region of the forcing grid.
!
!  CONTACT      : jan.polcher@lmd.jussieu.fr
!
!  LICENCE      : IPSL (2016)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!! 
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_driver/forcing_tools.f90 $
!! $Date: 2026-02-19 13:45:26 +0100 (jeu. 19 févr. 2026) $
!! $Revision: 9382 $
!_ ================================================================================================================================
!!
MODULE forcing_tools
  !
  USE defprec
  USE netcdf
  !
  USE ioipsl
  USE constantes
  USE solar
  !
  USE mod_orchidee_para
  USE forcingdaily_tools
  USE grid
  !
  IMPLICIT NONE
  ! 
  PRIVATE
  PUBLIC :: forcing_open, forcing_close, forcing_printdate, forcing_getvalues, forcing_printpoint,&
       &    forcing_getglogrid, forcing_givegridsize, forcing_givegrid, forcing_zoomgrid, forcing_integration_time, &
       &    forcing_meanprecip
      
  !
  !
  !
  INTERFACE forcing_reindex
     MODULE PROCEDURE forcing_reindex3d, forcing_reindex2dt, forcing_reindex2d, forcing_reindex1d, &
          &           forcing_reindex2to1, forcing_reindex1to2  
  END INTERFACE forcing_reindex
  !
  INTERFACE forcing_printpoint
     MODULE PROCEDURE forcing_printpoint_forgrid, forcing_printpoint_forgrid2d, forcing_printpoint_gen
  END INTERFACE forcing_printpoint
  !
  INTERFACE forcing_getvalues
     MODULE PROCEDURE forcing_getvalues1d, forcing_getvalues2d
  END INTERFACE forcing_getvalues
  !
  ! This PARAMETER essentially manages the memory usage of the module as it
  ! determines how much of the forcing will be uploaded from the netCDF file into
  ! memory.
  !
  INTEGER(i_std), SAVE :: slab_size_max
  !
  ! Time variables, all in Julian days
  !
  INTEGER(i_std), PARAMETER :: nbtmethods=4   ! number of cellmethods for forcing
  INTEGER(i_std), SAVE :: nbtax
  INTEGER(i_std), SAVE :: nb_forcing_steps
  REAL(r_std), SAVE :: global_start_date, global_end_date, forcing_tstep_ave
  REAL(r_std), SAVE :: dt_sechiba_keep
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)     :: time
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:)   :: time_bounds
  CHARACTER(LEN=20), SAVE, ALLOCATABLE, DIMENSION(:) :: time_axename, time_cellmethod
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)    :: time_sourcefile
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: time_id
  LOGICAL, SAVE :: end_of_file
  !
  ! Forcing file information
  !
  INTEGER(i_std), SAVE                                :: nb_forcefile=0
  CHARACTER(LEN=100), SAVE, ALLOCATABLE, DIMENSION(:) :: forfilename
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:), public  :: force_id
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: id_unlim
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: nb_atts, ndims, nvars
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)        :: convtosec
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: nbtime_perfile
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)      :: date0_file
  REAL(r_std), SAVE                                   :: startdate, forcingstartdate
  !
  ! Descrition of global grid
  !
  INTEGER(i_std), SAVE :: iim_glo, jjm_glo, nbpoint_glo, nbland_glo, nbseg
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: lon_glo, lat_glo
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:):: mask_glo
  INTEGER(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:):: area_glo
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)  :: lindex_glo
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: contfrac_glo
  LOGICAL, SAVE                                    :: compressed
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:):: corners_glo
  !
  ! Descritpion of zoomed grid
  !
  LOGICAL, SAVE :: zoom_forcing = .FALSE.
  INTEGER(i_std), SAVE :: iim_loc, jjm_loc, nbpoint_loc, nbland_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: lon_loc, lat_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)   :: lindex_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: area_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: contfrac_loc
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:):: corners_loc
  ! Number of land points per proc
  INTEGER(i_std), SAVE :: nbpoint_proc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:) :: glolindex_proc
  LOGICAL, SAVE :: landonly
  !-
  !- Heigh controls and data 
  !- 
  LOGICAL, SAVE                            :: zfixed, zsigma, zhybrid, zlevels, zheight 
  LOGICAL, SAVE                            :: zsamelev_uv 
  REAL, SAVE                               :: zlev_fixed, zlevuv_fixed 
  REAL, SAVE                               :: zhybrid_a, zhybrid_b 
  REAL, SAVE                               :: zhybriduv_a, zhybriduv_b
  LOGICAL, SAVE                            :: lwdown_cons
  !
  ! Forcing variables to be read and stored 
  !
  ! At 3000 we can fit in the slab an entire year of 3 hourly forcing.
  INTEGER(i_std), SAVE :: slab_size=-1
  INTEGER(i_std), SAVE :: current_offset=1
  INTEGER(i_std), SAVE :: position_slab(2)
  CHARACTER(LEN=20), SAVE :: calendar
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: tair_slab, qair_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: tairmax_slab, tairmin_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: hurs_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_tair, time_qair
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_tair, timebnd_qair
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: rainf_slab, snowf_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_precip
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_precip
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: preciptime_slab             !! Variable needed to keep track of how much rainfall was already distributed
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: swdown_slab, lwdown_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_swdown, time_lwdown
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_swdown, timebnd_lwdown
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: u_slab, v_slab, ps_slab
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_u, time_v, time_ps
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: timebnd_u, timebnd_v, timebnd_ps
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: ztq_slab, zuv_slab
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:)   :: reindex_glo, reindex_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: reindex2d_loc
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: origind
  !
  INTEGER(i_std), SAVE                              :: ncdfstart, ncdfcount
  !
CONTAINS
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_integration_time
!!
!>\BRIEF   Computes the interval over which the simulation should be carried out   
!!
!! DESCRIPTION:	 This routing will get the following parameters from the run.def : 'START_DATE', 'END_DATE' and 'DT_SECHIBA'.
!!               It allows to define the integration time of ORCHIDEE and later it will be used to verify that we have
!!               the needed data in the forcing files to perform this simulation.
!!
!! \n
!_ ==============================================================================================================================
!!
  SUBROUTINE forcing_integration_time(date_start, dt, nbdt)
    !
    !
    ! This subroutine gets the start date of the simulation, the time step and the number 
    ! of time steps we need to do until the end of the simulations.
    !
    !
    !
    REAL(r_std), INTENT(out)                     :: date_start     !! The date at which the simulation starts
    REAL(r_std), INTENT(out)                     :: dt             !! Time step length in seconds
    INTEGER(i_std), INTENT(out)                  :: nbdt           !! Number of timesteps to be executed
    !
    ! Local
    !
    CHARACTER(LEN=20) :: str_sdate(2), str_edate(2), tmpstr
    INTEGER(i_std) :: s_year, s_month, s_day, e_year, e_month, e_day
    INTEGER(i_std) :: seci, hours, minutes
    REAL(r_std) :: s_sec, e_sec, dateend, diff_sec, date_end
    INTEGER(i_std) :: i, ic
    CHARACTER(LEN=20) :: str_cyclic_start(2), str_cyclic_end(2)
    INTEGER(i_std) :: c_year_start, c_month_start, c_day_start, c_year_end, c_month_end, c_day_end
    REAL(r_std) :: c_sec_start, c_sec_end
    !
    !Config Key  = START_DATE
    !Config Desc = Date at which the simulation starts
    !Config Def  = NONE
    !Config If   = orchideedriver
    !Config Help = The format is the same as in the CF convention : 1999-09-13 12:0:0
    !Config Units= [yyyy-mm-dd hh:0:0]
    str_sdate = " "
    CALL getin('START_DATE',str_sdate)
    !
    !Config Key  = CYCLIC_STARTDATE
    !Config Desc = Date at which the cyclic year is started
    !Config Def  = NONE
    !Config If   = orchideedriver
    !Config Help = The format is the same as in the CF convention : 1999-09-13 12:0:0
    !Config Units= [yyyy-mm-dd hh:0:0]
    str_cyclic_start = "NONE"
    CALL getin('CYCLIC_STARTDATE',str_cyclic_start)
    
    !
    !Config Key  = CYCLIC_ENDDATE
    !Config Desc = Date at which the cyclic year is ended
    !Config Def  = NONE
    !Config If   = orchideedriver
    !Config Help = The format is the same as in the CF convention : 1999-09-13 12:0:0
    !Config Units= [yyyy-mm-dd hh:0:0]
    str_cyclic_end = "NONE"
    CALL getin('CYCLIC_ENDDATE',str_cyclic_end)
    !
    ! the start date of simulation 
    IF ( (INDEX(str_sdate(1),"-") .NE. INDEX(str_sdate(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_sdate(2),":") .NE. INDEX(str_sdate(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_sdate(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_sdate(1) = tmpstr
          tmpstr = str_sdate(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_sdate(2) = tmpstr
       ENDDO
       READ (str_sdate(1),*) s_year, s_month, s_day
       READ (str_sdate(2),*) hours, minutes, seci
       s_sec = hours*3600. + minutes*60. + seci
    ELSE
       CALL ipslerr(3, "forcing_integration_time", "START_DATE incorrectly specified in run.def", str_sdate(1), str_sdate(2))
    ENDIF
    !---------------------------------
    ! cyclic start date
    IF (str_cyclic_start(1) .EQ. "NONE") THEN
       ! CYCLIC_STARTDATE is not set in run.def
       ! Set same values as START_DATE
       c_year_start=s_year
       c_month_start=s_month
       c_day_start=s_day
       c_sec_start=s_sec

    ELSE IF ( (INDEX(str_cyclic_start(1),"-") .NE. INDEX(str_cyclic_start(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_cyclic_start(2),":") .NE. INDEX(str_cyclic_start(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_cyclic_start(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_cyclic_start(1) = tmpstr
          tmpstr = str_cyclic_start(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_cyclic_start(2) = tmpstr
       ENDDO
       READ (str_cyclic_start(1),*) c_year_start, c_month_start, c_day_start
       READ (str_cyclic_start(2),*) hours, minutes, seci
       c_sec_start = hours*3600. + minutes*60. + seci
    ELSE IF ( len_trim(str_cyclic_start(1)) .NE. 0 ) THEN
       CALL ipslerr(3, "forcing_integration_time", "CYCLIC_STARTDATE incorrectly specified in run.def", str_cyclic_start(1), str_cyclic_start(2))
    ENDIF
    ! if s_year not the same as c_year, use cyear to compute date_start
    IF ( ( s_year .NE. c_year_start) .AND. (len_trim(str_cyclic_start(1)) .NE. 0)) THEN
       CALL ymds2ju (c_year_start, c_month_start, c_day_start, c_sec_start, date_start)
    ELSE
       CALL ymds2ju (s_year, s_month, s_day, s_sec, date_start)
    ENDIF
    CALL forcing_printdate(date_start, "This is after reading the start date")

    !
    !Config Key  = END_DATE
    !Config Desc = Date at which the simulation ends
    !Config Def  = NONE
    !Config If   = orchideedriver
    !Config Help = The format is the same as in the CF convention : 1999-09-13 12:0:0
    !Config Units= [yyyy-mm-dd hh:0:0]
    str_edate = " "
    CALL getin('END_DATE',str_edate)
    !
    IF ( (INDEX(str_edate(1),"-") .NE. INDEX(str_edate(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_edate(2),":") .NE. INDEX(str_edate(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_edate(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_edate(1) = tmpstr
          tmpstr = str_edate(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_edate(2) = tmpstr
       ENDDO
       READ (str_edate(1),*) e_year, e_month, e_day
       READ (str_edate(2),*) hours, minutes, seci
       e_sec = hours*3600. + minutes*60. + seci
    ELSE
       CALL ipslerr(3, "forcing_integration_time", "END_DATE incorrectly specified in run.def", str_edate(1), str_edate(2))
    ENDIF

    !---------------------------------
    ! for cyclic end date
    IF (str_cyclic_end(1) .EQ. "NONE") THEN
       ! CYCLIC_ENDDATE is not set in run.def
       ! Set same values as END_DATE
       c_year_end=e_year
       c_month_end=e_month
       c_day_end=e_day
       c_sec_end=e_sec
    ELSE IF ( (INDEX(str_cyclic_end(1),"-") .NE. INDEX(str_cyclic_end(1),"-", .TRUE.)) .AND. &
         &  (INDEX(str_cyclic_end(2),":") .NE. INDEX(str_cyclic_end(2),":", .TRUE.)) ) THEN
       DO i=1,2
          tmpstr = str_cyclic_end(1)
          ic = INDEX(tmpstr,"-")
          tmpstr(ic:ic) = " "
          str_cyclic_end(1) = tmpstr
          tmpstr = str_cyclic_end(2)
          ic = INDEX(tmpstr,":")
          tmpstr(ic:ic) = " "
          str_cyclic_end(2) = tmpstr
       ENDDO
       READ (str_cyclic_end(1),*) c_year_end, c_month_end, c_day_end
       READ (str_cyclic_end(2),*) hours, minutes, seci
       c_sec_end = hours*3600. + minutes*60. + seci
    ELSE IF ( len_trim(str_cyclic_end(1)) .NE. 0 ) THEN
       CALL ipslerr(3, "forcing_integration_time", "CYCLIC_ENDDATE incorrectly specified in run.def", str_cyclic_end(1), str_cyclic_end(2))
    ENDif

    ! if e_year not the same as c_year_end, use cyear_end to compute date_end
    IF (( e_year .NE. c_year_end)  .AND. (len_trim(str_cyclic_end(1)) .NE. 0) )THEN
       CALL ymds2ju (c_year_end, c_month_end, c_day_end, c_sec_end, date_end)
    ELSE
       CALL ymds2ju (e_year, e_month, e_day, e_sec, date_end)
    ENDIF
    
    !
    IF (( s_year .NE. c_year_start) .AND. (len_trim(str_cyclic_start(1)) .NE. 0) )then
       CALL time_diff (c_year_start,c_month_start,c_day_start,c_sec_start,c_year_end,c_month_end,c_day_end,c_sec_end,diff_sec)
    ELSE
       CALL time_diff (s_year,s_month,s_day,s_sec,e_year,e_month,e_day,e_sec,diff_sec)
    ENDIF

    !
    !Config Key  = DT_SECHIBA
    !Config Desc = Time step length in seconds for sechiba component
    !Config Def  = 1800
    !Config If   = orchideedriver
    !Config Help = 
    !Config Units= [seconds]
    dt = 1800
    CALL getin('DT_SECHIBA', dt)
    dt_sechiba_keep = dt
    !
    nbdt = NINT(diff_sec/dt)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Read the configuration options for the time interpolations.
    !
    !Config Key   = LWDOWN_CONS
    !Config Desc  = Conserve the longwave downward radiation of the forcing
    !Config Def   = n
    !Config If    = orchideedriver
    !Config Help  = This flag allows to conserve the downward longwave radiation
    !               provided in the forcing. It will do this by taking the closest
    !               neighbour in time from the forcing. This assumes that the forcing
    !               contains average fluxes. The default setting (LWDOWN_CONS=n) will
    !               force the model to perform a linear interpolation of the fluxes.
    !Config Units = [y/n]
    !-
    lwdown_cons = .FALSE.
    CALL getin('LWDOWN_CONS', lwdown_cons)
    !
  END SUBROUTINE forcing_integration_time
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_open
!!
!>\BRIEF      Opens the forcing files and extract the main information.
!!
!! DESCRIPTION:	 This routine opens all the forcing files provided in the list and verifies that the grid corresponds
!!               to the coordinates provided (and which was obtained by the model from glogrid.f90.). It then zooms
!!               into the forcing as requested by the user, extracts the vertical coordinates and final reads the time axis.
!!               Some basic consistency checks are performed as for instance ensuring the that all the forcing data is available
!!               to simulate the desired period.
!!               All that information is also broadcasted to all processors.
!!               Actual forcing data is not read at this stage.
!!
!! \n
!_ ==============================================================================================================================
!
  SUBROUTINE forcing_open(filenames_in, iim, jjm, lon, lat, nbpoint_in, drvzoom_lon, drvzoom_lat, &
       &                  kindex, nbindex_perproc, wunit, model_guess, zoom_index_init, landonly_arg)
    !
    ! Opens the forcing file and reads some key information and stores them in the shared variables of the
    ! module.
    !
    ! Lon, lat should come from the grid file we read before. This will give indication of the grid
    ! file is consistant with the forcing file and if we need to zoom into the forcing file.
    !
    ! Time interval of the simulation is also determined.
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in)          :: filenames_in(:)
    INTEGER(i_std), INTENT(in)            :: iim, jjm, nbpoint_in
    REAL(r_std), INTENT(in)               :: lon(iim,jjm), lat(iim,jjm)
    REAL(r_std), DIMENSION(2), INTENT(in) :: drvzoom_lon, drvzoom_lat
    INTEGER(i_std), INTENT(in)            :: kindex(nbpoint_in)
    INTEGER(i_std), INTENT(in)            :: nbindex_perproc
    INTEGER(i_std), OPTIONAL, INTENT(in)  :: wunit
    LOGICAL, OPTIONAL, INTENT(in)         :: landonly_arg
    CHARACTER(LEN=*), INTENT(in)          :: model_guess
    INTEGER, DIMENSION(2), INTENT(in)     :: zoom_index_init
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iim_tmp, jjm_tmp, nbpoint_tmp, nb_files    
    INTEGER(i_std) :: iv, it
    INTEGER(i_std) :: inl, ii, jj, ik
    INTEGER(i_std) :: land_id
    REAL(r_std)    :: dt
    INTEGER(i_std) :: nbdt
    LOGICAL        :: if_dump
    INTEGER, DIMENSION(2) :: iim_jjm_grid
    !
    ! Check optional arguments
    !
    ! The default behaviour is to provide only land points to the calling program.
    ! But for forcing ocean model there is also the option to pass on land and ocean values.
    ! When the grid is initialized landonly_tmp=.FALSE. has to be set to obtian this behaviour.
    !
    IF ( PRESENT(landonly_arg) ) THEN
       landonly=landonly_arg
    ELSE IF ( model_guess .EQ. 'unstruct' ) THEN
       landonly = .FALSE.
    ELSE
       landonly = .TRUE.
    ENDIF

    ! Predefine if_dump to false, that is not to dump forcing data into a temporary file in RUN_DIR
    ! this temporary file would be useful for verification purpose but doesn't work with unstructured grid.
    if_dump = .FALSE.
    
    !Config Key  = FORCING_MEMORY
    !Config Desc = Number of time steps of the forcing we will store in memory.
    !Config Def  = 80
    !Config If   = orchideedriver
    !Config Help = To reduce and optimise disk accesses more of the forcing can be loaded into
    !Config        memory. With this parameter the amount of memory can be adjusted. Be carefull
    !Config        as if you use too much memory the system will cick you out or slow down the
    !Config        execution of your program. 
    !Config Units= [-]
    slab_size_max = 80
    CALL getin('FORCING_MEMORY', slab_size_max)
    !
    ! How many files do we have to open ?
    !
    !
    ! All the meta information from the forcing file is ojnly needed on the root processor.
    !
    IF ( is_root_prc ) THEN
       !
       CALL forcing_filenamecheck(filenames_in, nb_files)
       IF ( PRESENT(wunit) ) THEN
          DO it=1,nb_files
             WRITE(wunit,*) "Files to be used for forcing the simulation :", it, TRIM(forfilename(it))
          ENDDO
       ENDIF
       !
       ! 0.0 Check if variables are allocated to the right size on root_proc
       !
       IF (nb_files > nb_forcefile) THEN
          IF ( ALLOCATED(force_id) ) DEALLOCATE(force_id)
          ALLOCATE(force_id(nb_files))
          IF ( ALLOCATED(id_unlim) )  DEALLOCATE(id_unlim)
          ALLOCATE(id_unlim(nb_files))
          IF ( ALLOCATED(nb_atts) ) DEALLOCATE(nb_atts)
          ALLOCATE(nb_atts(nb_files))
          IF ( ALLOCATED(ndims) ) DEALLOCATE(ndims)
          ALLOCATE(ndims(nb_files))
          IF ( ALLOCATED(nvars) ) DEALLOCATE(nvars)
          ALLOCATE( nvars(nb_files))
          IF ( ALLOCATED(nbtime_perfile) ) DEALLOCATE(nbtime_perfile)
          ALLOCATE(nbtime_perfile(nb_files))
          IF ( ALLOCATED(convtosec) ) DEALLOCATE(convtosec)
          ALLOCATE(convtosec(nb_files))
       ENDIF
       nb_forcefile = nb_files
       !
       ! Get the global grid size from the forcing file. The output is in temporary variables as in this
       ! module the values are shared.
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting global grid from ",  nb_forcefile, "files."
          CALL FLUSH(wunit)
       ENDIF
       !
       ! for unstructured grids, if_dump must be false
       IF ( model_guess .EQ. 'unstruct'  .AND. if_dump) THEN
          CALL ipslerr (3,'forcing_open',"for unstructured grids, if_dump must be false."," "," ")
       ENDIF
       !
       CALL forcing_getglogrid(nb_forcefile, forfilename, iim_tmp, jjm_tmp, nbpoint_tmp, if_dump, landonly)
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting the full grid", nbpoint_tmp
          CALL FLUSH(wunit)
       ENDIF
       !
       iim_jjm_grid = (/iim,jjm/)
       !CALL forcing_zoomgrid(drvzoom_lon, drvzoom_lat, forfilename(1), model_guess, if_dump, iim, jjm, zoom_index_init)
       CALL forcing_zoomgrid(drvzoom_lon, drvzoom_lat, forfilename(1), model_guess, if_dump, iim_jjm_grid, zoom_index_init)
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Out of the zoomed grid operation"
          CALL FLUSH(wunit)
       ENDIF
       !
       ! Verification that the grid sizes coming from the calling program are consistant with what we get 
       ! from the forcing file. iim, jjm: from grid file. iim_loc from forcing file.
       !
       IF ( (iim_loc .NE. iim) .OR. (jjm_loc .NE. jjm) ) THEN
          CALL ipslerr (3,'forcing_open',"At least one of the dimensions of the grid obtained from the",&
               &        "grid file is different from the one in the forcing file.",&
               &        "Run driver2oasis -init to generate a new grid file.")
       ENDIF
       ! Special treatment for the number of land point, as we could have a case where the forcing
       ! file does not include the land/sea mask.
       !
       IF ( nbpoint_loc .NE. nbpoint_in ) THEN
          ! We trust the number of land points obtained from the gridfile. It has the land/sea mask.
          nbpoint_loc = nbpoint_in
       ENDIF
       
       !
       ! Treat the time dimension now :
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting forcing time"
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_time(nb_forcefile, forfilename)
       !
       ! Now that we know how much time steps are in the forcing we can set some realistic slab_size
       !
       slab_size=MIN(nb_forcing_steps, slab_size_max)
       !
       !
       ! Get the vertical information from the file
       !
       CALL forcing_vertical(force_id(1))
       !
       !
       IF ( PRESENT(wunit) ) THEN
          WRITE(wunit,*) "Getting integration time"
          CALL FLUSH(wunit)
       ENDIF
       CALL forcing_integration_time(startdate, dt, nbdt)

       ! Test that the time interval requested by the user correspond to the time available in the 
       ! forcing file.
       !
       IF ( startdate < time_bounds(1,1,1) .OR. startdate > time_bounds(nb_forcing_steps,1,2) ) THEN
          CALL forcing_printdate(startdate, "--> Start Date in forcing_open")
          CALL forcing_printdate(time_bounds(1,1,1), "--> Outer bound of forcing file.")
          CALL forcing_printdate(time_bounds(nb_forcing_steps,1,2), "--> Last date to be simulated.")
          CALL ipslerr (3,'forcing_open', 'Start time requested by the user is outside of the time interval',&
               & "covered by the forcing file.","Please verify the configuration in the run.def file.")
          
       ENDIF
       !
       IF ( startdate+(dt/one_day)*nbdt > time_bounds(nb_forcing_steps,1,2) .OR. &
            & startdate+(dt/one_day)*nbdt < time_bounds(1,1,1)) THEN
          CALL forcing_printdate(time_bounds(nb_forcing_steps,1,2), "Outer bound of forcing file.")
          CALL forcing_printdate(startdate+(dt/one_day)*nbdt, "Last date to be simulated.")
          WRITE(*,*) "ERROR : Final date of forcing needed is : ", startdate+(dt/one_day)*nbdt
          WRITE(*,*) "ERROR : The outer bound of the last forcing time step is :", time_bounds(nb_forcing_steps,1,2)
          CALL ipslerr (3,'forcing_open', 'End time requested by the user is outside of the time interval',&
               & "covered by the forcing file.","Please verify the configuration in the run.def file.")
       ENDIF
       !
    ENDIF
    !
    ! Broadcast the local grid (i.e. the one resulting from the zoom) to all processors
    !
    CALL bcast(iim_loc)
    CALL bcast(jjm_loc)
    CALL bcast(nbpoint_loc)
    CALL bcast(nbland_loc)
    ! Time variables needed by all procs
    CALL bcast(slab_size)
    CALL bcast(startdate)
    CALL bcast(forcingstartdate)
    CALL bcast(forcing_tstep_ave)
    CALL bcast(nbseg)
    !
    ! Number of points per processor
    !
    IF ( landonly .OR. model_guess .EQ. 'unstruct' ) THEN
       nbpoint_proc = nbindex_perproc
    ELSE
       nbpoint_proc = nbpoint_glo
    ENDIF
    !
    ! On the slave processes we need to allocate the memory for the data on root_prc to be bcast
    ! On the root_proc these allocations were done with CALL forcing_zoomgrid
    !
    ALLOCATE(glolindex_proc(nbpoint_proc))
    IF ( .NOT. is_root_prc ) THEN
       ALLOCATE(lon_loc(iim_loc,jjm_loc))
       ALLOCATE(lat_loc(iim_loc,jjm_loc))
       ALLOCATE(lindex_loc(nbpoint_loc)) 
       ALLOCATE(mask_loc(iim_loc,jjm_loc))
       ALLOCATE(area_loc(iim_loc,jjm_loc))
       ALLOCATE(contfrac_loc(nbpoint_loc))
       ALLOCATE(corners_loc(iim_loc,jjm_loc,nbseg,2))
    ENDIF
    !
    ! Keep on each processor the index of each land point on the *_loc grid
    !
    IF ( landonly .OR. model_guess .EQ. 'unstruct') THEN
       CALL scatter(kindex, glolindex_proc)
    ELSE
       !
       ! Build a simple indexing list as the one for land cannot be used.
       !
       ik=0
       DO jj=1,jjm_loc
          DO ii=1,iim_loc
             ik=ik+1
             glolindex_proc(ik) = ik
          ENDDO
       ENDDO
    ENDIF
    !
    CALL bcast(lon_loc)
    CALL bcast(lat_loc)
    CALL bcast(lindex_loc)
    CALL bcast(mask_loc)
    CALL bcast(area_loc)
    CALL bcast(contfrac_loc)
    CALL bcast(corners_loc)
    !
  END SUBROUTINE forcing_open
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_getvalues1d
!!
!>\BRIEF   Gets the forcing data for a time interval of simulation (usually 30min interval). 1d means 1D in space.  
!!
!! DESCRIPTION:	The routine will get the forcing valid for the time interval provided by the caller.
!!              First it will check that the data is already in memory for that time interval. If not
!!              it will first read the data from the netCDF file.
!!              Then the forcing date will be interpolated to the requested time interval.
!!              The code calls linear interpolation for most variables except for SWdown and precipitation.
!!              These temporal interpolations can be improved later.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_getvalues1d(time_int, dt, zlev_tq, zlev_uv, tair, qair, rainf, snowf, &
       &                       swdown, lwdown, solarang, u, v, ps, model_guess, is_spinupacc)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)  :: dt                                     !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    CHARACTER(LEN=8), OPTIONAL, INTENT(in) :: model_guess
    LOGICAL, OPTIONAL, INTENT(in) :: is_spinupacc
    
    REAL(r_std), INTENT(out) :: zlev_tq(:), zlev_uv(:)
    REAL(r_std), INTENT(out) :: tair(:), qair(:), rainf(:), snowf(:)
    REAL(r_std), INTENT(out) :: swdown(:), lwdown(:), solarang(:)
    REAL(r_std), INTENT(out) :: u(:), v(:), ps(:)
    
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i
    CHARACTER(LEN=8) :: tmp_modelguess
    LOGICAL :: is_tmp_spinup
    !
    IF ( PRESENT(model_guess)) THEN
       tmp_modelguess = model_guess 
    ELSE
      CALL ipslerr (3,'forcing_getvalues1d', 'model_guess is not given',&
           & "","")
    ENDIF
       
    IF ( PRESENT(is_spinupacc)) THEN
       is_tmp_spinup = is_spinupacc
    ELSE
       is_tmp_spinup = .FALSE.
    ENDIF
    !
    ! Test that we have the time interval within our slab of data else we need to update it.
    ! Att : the tests are done here on time_tair as an exemple. This might need to have to be generalized.
    !
    !
    ! First case the time axis of the variable are not even yet allocated !
    !
    IF ( .NOT. ALLOCATED(time_tair) ) THEN
       ! read in the slab of forcing file which correspond to the time_int interval, and the slab size
       ! these data are then distributed to different processors, with variable names *_slab)
       ! in fact, slabsize could be much smaller than 80, because the time interval is only 30 min.
       ! But more tests can be done to see the impact of different slabsize on the computing speed
       CALL forcing_readslab(time_int, model_guess,is_tmp_spinup)
       CALL forcing_printdate(timebnd_tair(1,1), "Start of time slab just read", numout)
       CALL forcing_printdate(time_tair(1), "Time of first temperature value", numout)
       CALL forcing_printdate(timebnd_tair(slab_size,2), "End of time slab just read", numout)
    ELSE
       ! If we have time axis (for TAIR here) we test that it is long enough in time to allow for an interpolation.
       !
       IF ( time_int(2)+forcing_tstep_ave/one_day > time_tair(slab_size) .AND. (.NOT. end_of_file) ) THEN
          CALL forcing_readslab(time_int, model_guess,is_tmp_spinup)
          CALL forcing_printdate(timebnd_tair(1,1), "Start of time slab just read", numout)
          CALL forcing_printdate(time_tair(1), "Time of first temperature value", numout)
          CALL forcing_printdate(timebnd_tair(slab_size,2), "End of time slab just read", numout)
       ENDIF
    ENDIF
    !
    ! temporal interpolation of forcing data to the simulation time step
    ! in unstructured dynamico grids, the forcing_tstep_ave is 3 hours
    IF ( forcing_tstep_ave <= one_day/3.0) THEN
       !
       ! Interpolate the dynamical variables to the time step at which the driver is for the moment.
       !
       CALL forcing_interpol(time_int, dt, time_u, u_slab, u)
       CALL forcing_interpol(time_int, dt, time_v, v_slab, v)
       CALL forcing_interpol(time_int, dt, time_ps, ps_slab, ps)
       !
       ! Compute the height of the first level (a routine will be needed for that !)
       ! ATT : we assume that the time axis for the height of the scalar variable is the one of TAIR
       ! and for the height of wind is the same as U.
       
       CALL forcing_interpol(time_int, dt, time_tair, ztq_slab, zlev_tq)
       CALL forcing_interpol(time_int, dt, time_u, zuv_slab, zlev_uv)
       !
       ! Interpolate the state variables of the lower atmospheric level
       !
       CALL forcing_interpol(time_int, dt, time_tair, tair_slab, tair)
       
       CALL forcing_interpol(time_int, dt, time_qair, qair_slab, qair)
       !
       ! Spread the precipitation as requested by the user
       !
       CALL forcing_spreadprec(time_int, dt, timebnd_precip, time_precip, rainf, snowf)
       !
       ! Deal with the interpolate of the radiative fluxes.
       CALL forcing_solarint(time_int, dt, timebnd_swdown, time_swdown, iim_loc, &
            jjm_loc, lon_loc, lat_loc, tmp_modelguess, swdown, solarang)
       !
       ! We have the option here to conserve LWdown by taking the closest point in the forcing.
       ! So no interpolation is done.
       !
       IF ( lwdown_cons ) THEN
          CALL forcing_closest(time_int, dt, time_lwdown, lwdown_slab, lwdown)
       ELSE
          CALL forcing_interpol(time_int, dt, time_lwdown, lwdown_slab, lwdown)
       ENDIF
       !
    ELSE IF (forcing_tstep_ave == one_day) THEN
       !
       ! If the forcing is daily we move to the module designed for these interpolations
       ! that is, to generate the diurnal cycle according to dt and based on the Tmax and Tmin
       !
       CALL forcingdaily_gensubd(time_int, dt, iim_loc, jjm_loc, lon_loc, lat_loc, glolindex_proc, &
            &                    nbpoint_proc, slab_size, time_tair, ztq_slab, zuv_slab, tair_slab, &
            &                    tairmin_slab, tairmax_slab, qair_slab, hurs_slab, rainf_slab, &
            &                    snowf_slab, swdown_slab, lwdown_slab, u_slab, v_slab, ps_slab)
       !
       ! Extracts the forcing values needed by ORCHIDEE from the re-generated diurnal cycles. Based
       ! on the description of the subroutine these are the values that will be used to help in
       ! the regeneration of the diurnal cycle of the next time step.
       !
       CALL forcingdaily_getvalues(time_int, dt, zlev_tq, zlev_uv, tair, qair, rainf, snowf, &
            &                      swdown, lwdown, solarang, u, v, ps)
       !
    ELSE
       !
       ! Catch any forcing files not adapted to the interpolations available.
       !
       WRITE(numout,*) "#### Forcing time step is too long for a the foreseen interpolations"
       WRITE(numout,*) "#### forcing_tstep_ave =", forcing_tstep_ave
       CALL ipslerr (3,'forcing_getvalues1d', 'Forcing time step incompatible with implemented interpolations.',&
            & "","")
    ENDIF
    !
  END SUBROUTINE forcing_getvalues1d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_getvalues2d
!!
!>\BRIEF   Gets the forcing data in 2D field for a time interval. 2D means 2D in space. 
!!
!! DESCRIPTION:	The routine will get the forcing valid for the time interval provided by the caller.
!!              First it will check that the data is already in memory for that time interval. If not
!!              it will first read the data from the netCDF file.
!!              Then the forcing date will be interpolated to the requested time interval.
!!              The code calls linear interpolation for most variables except for SWdown and precipitation.
!!              These temporal interpolations can be improved later.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_getvalues2d(time_int, dt, zlev_tq, zlev_uv, tair, qair, rainf, snowf, &
       &                       swdown, lwdown, solarang, u, v, ps, model_guess, is_spinupacc)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)  :: dt                                     !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    CHARACTER(LEN=8), OPTIONAL, INTENT(in) :: model_guess
    LOGICAL, OPTIONAL, INTENT(in) :: is_spinupacc
    
    REAL(r_std), INTENT(out) :: zlev_tq(:,:), zlev_uv(:,:)
    REAL(r_std), INTENT(out) :: tair(:,:), qair(:,:), rainf(:,:), snowf(:,:)
    REAL(r_std), INTENT(out) :: swdown(:,:), lwdown(:,:), solarang(:,:)
    REAL(r_std), INTENT(out) :: u(:,:), v(:,:), ps(:,:)
    !
    REAL(r_std) :: zzlev_tq(nbpoint_loc), zzlev_uv(nbpoint_loc)
    REAL(r_std) :: ztair(nbpoint_loc), zqair(nbpoint_loc), zrainf(nbpoint_loc), zsnowf(nbpoint_loc)
    REAL(r_std) :: zswdown(nbpoint_loc), zlwdown(nbpoint_loc), zsolarang(nbpoint_loc)
    REAL(r_std) :: zu(nbpoint_loc), zv(nbpoint_loc), zps(nbpoint_loc)
    INTEGER(i_std) :: i, j, k
    !
    CALL forcing_getvalues(time_int, dt, zzlev_tq, zzlev_uv, ztair, zqair, zrainf, zsnowf, zswdown, zlwdown, zsolarang, zu, zv, zps)
    !
    k = 0
    DO j=1,jjm_loc
       DO i=1,iim_loc
          k = k + 1
          zlev_tq(i,j) = zzlev_tq(k)
          zlev_uv(i,j) = zzlev_uv(k)
          tair(i,j) = ztair(k)
          qair(i,j) = zqair(k)
          rainf(i,j) = zrainf(k)
          snowf(i,j) = zsnowf(k)
          swdown(i,j) = zswdown(k)
          lwdown(i,j) = zlwdown(k)
          solarang(i,j) = zsolarang(k)
          u(i,j) = zu(k)
          v(i,j) = zv(k)
          ps(i,j) = zps(k)
       ENDDO
    ENDDO
    !
    !
  END SUBROUTINE forcing_getvalues2d
    
!!  =============================================================================================================================
!! SUBROUTINE: forcing_closest
!!
!>\BRIEF   This routine does not interpolate and simply uses the closes value in time. It is useful for preserving
!!         variables which are averaged in the forcing file.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_closest(time_int_in, dt, time_central_in, var_slab, var)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int_in(2)
    REAL(r_std), INTENT(in)  :: dt
    REAL(r_std), INTENT(in)  :: time_central_in(:)
    REAL(r_std), INTENT(in)  :: var_slab(:,:)
    REAL(r_std), INTENT(out) :: var(:)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: slabind_a, slabind_b, imin(1), i
    REAL(r_std) :: time_int(2), time_central(slab_size_max)
    REAL(r_std) :: mid_int, wa, wb, wt, wab, wae, tmp_mid_int
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
    !
    ! Shift the input dates in order to gain in precision for the calculations 
    !
    IF ( .NOT. ALLOCATED(mask) ) THEN
       ALLOCATE(mask(slab_size_max))
       mask(:) = .FALSE.
    ENDIF
    !
    time_int(:) = time_int_in(:)-INT(forcingstartdate)
    time_central(1:slab_size) = time_central_in(1:slab_size)-INT(forcingstartdate)
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    mask(1:slab_size) = .TRUE.
    !
    ! Select the forcing interval for which the center date is the closest to the time of 
    ! the model.
    !
    mid_int = time_int(1) + (dt/2.0)/one_day
    imin = MINLOC( ABS(time_central(1:slab_size) - mid_int), mask(1:slab_size) )
    !
    ! Verify that this is a possible date
    !
    IF ( imin(1) > 0 .AND. imin(1) <= slab_size ) THEN
       !
       slabind_a = imin(1)
       !
    ELSE
       CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
       CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
       CALL forcing_printdate(time_central_in(imin(1)), "===> Center of forcing time interval.")
       CALL ipslerr (3,'forcing_closest', 'The target time interval has no acceptable closest',&
            & "time in the forcing slab.","")
    ENDIF
    !
    ! Transfer the data from the sloest time of the forcing data slab.
    !
    DO i=1, nbpoint_proc
       !
       var(i) = var_slab(i,slabind_a)
       !
    ENDDO
    !
    !
  END SUBROUTINE forcing_closest
  
!!  =============================================================================================================================
!! SUBROUTINE: forcing_interpol
!!
!>\BRIEF   Perform linear interpolation for the time interval requested.
!!
!! DESCRIPTION:	  
!! The code gets an interval over which the model will integrate (time_int_in) but only uses the centre. It also gets 
!! the times representative of the forcing data for the variable at hand (time_central_in). Using this data we will  
!! determine which 2 forcing times will need to be used for the interpolation. Once this is established the weights 
!! are computed and used in order to interpolate the variable between the 2 times which bracket the model integration time. 
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_interpol(time_int_in, dt, time_central_in, var_slab, var)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int_in(2)             !! The time interval over which the forcing is needed by the model.
    REAL(r_std), INTENT(in)  :: dt                         !! Time step of the model
    REAL(r_std), INTENT(in)  :: time_central_in(:)         !! Representative time for the interval of validity of the forcing data
    REAL(r_std), INTENT(in)  :: var_slab(:,:)              !! The slab of forcing data read from the file.
    REAL(r_std), INTENT(out) :: var(:)                     !! Result of the time interpolation.
    !
    ! LOCAL
    !
    INTEGER(i_std) :: slabind_a, slabind_b, imin(1), i
    REAL(r_std) :: time_int(2), time_central(slab_size_max)
    REAL(r_std) :: mid_int, wa, wb, wt, wab, wae, tmp_mid_int
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    IF ( .NOT. ALLOCATED(mask) ) THEN
       ALLOCATE(mask(slab_size_max))
       mask(:) = .TRUE.
    ENDIF
    !
    ! Shift the input dates in order to gain in precision for the calculations 
    !
    time_int(:) = time_int_in(:)-INT(forcingstartdate)
    time_central(1:slab_size) = time_central_in(1:slab_size)-INT(forcingstartdate)
    !
    ! Select the type of interpolation to be done.
    !
    ! Compute the central time of the model integration time.
    !
    mid_int = time_int(1) + (dt/2.0)/one_day
    ! Locate that time on the time axis of the forcing.
    imin = MINLOC( ABS(time_central(1:slab_size) - mid_int), mask(1:slab_size) )
    !
    ! Determine which indices are to the left (slabind_a) and right (slabind_b) of the model time and will be used
    ! for the linear interpolation.
    !

    IF ( imin(1) > 1 .AND. imin(1) < slab_size ) THEN
       !
       ! Determine if the model time is to the left or right of the representative time
       ! of the forcing data. This allows to determine with which other position in the
       ! forcing data we need to interpolate.
       !
       IF ( mid_int < time_central(imin(1)) ) THEN
          slabind_a = imin(1) - 1
          slabind_b = imin(1)
       ELSE
          slabind_a = imin(1)
          slabind_b = imin(1) + 1
       ENDIF
       !
    ELSE IF ( imin(1) == 1 ) THEN
       !
       ! If we are at the first time step of the forcing data we need to take care as there is 
       ! no data earlier.
       !
       slabind_a = 1
       slabind_b = 2
       IF ( mid_int < time_central(slabind_a) ) THEN
          IF ( time_int(2) < time_central(slabind_a) ) THEN
             CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
             CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
             CALL forcing_printdate(time_central_in(slabind_a), "===> Center of forcing time interval.")
             CALL ipslerr (3,'forcing_interpol', 'The target time interval lies before the first date of the slab.',&
                  & "","")
          ELSE
             mid_int = time_central(slabind_a) 
          ENDIF
       ENDIF
    ELSE IF ( imin(1) == slab_size ) THEN
       !
       ! If we are at the end of the forcing data we need to pay attention as we have no data later in time.
       !
       slabind_a = slab_size - 1
       slabind_b = slab_size
       IF ( mid_int > time_central(slabind_b) ) THEN
          IF ( time_int(1) > time_central(slabind_b) ) THEN
             CALL forcing_printdate(time_int_in(1), "===> Start of target time interval.")
             CALL forcing_printdate(time_int_in(2), "===> End of target time interval.")
             CALL forcing_printdate(time_central_in(slabind_b), "===> Center of forcing time interval.")
             CALL ipslerr (3,'forcing_interpol', 'The target time interval lies after the last date of the slab.',&
                  & "","")
          ELSE
             mid_int = time_central(slabind_b) 
          ENDIF
       ENDIF
    ENDIF
    !
    ! Compute the weights between the values at slabind_a and slabind_b. As with the time
    ! representation we are at the limit of precision we use 2 days to compute the distance
    ! in time between the first value (slabind_a) and the middle of the target interval.
    !
    wab = time_int(1) - time_central(slabind_a) + (dt/2.0)/one_day
    wae = time_int(2) - time_central(slabind_a) - (dt/2.0)/one_day
    wa = (wab+wae)/2.0
    wb = time_central(slabind_b) - time_central(slabind_a)
    wt = wa/wb
    !
    ! Do the weighted average of all land points with the time indices and weights computed above.
    !
    DO i=1, nbpoint_proc
       var(i) = var_slab(i,slabind_a) + wt*(var_slab(i,slabind_b) - var_slab(i,slabind_a))
    ENDDO

  END SUBROUTINE forcing_interpol


!!  =============================================================================================================================
!! SUBROUTINE: forcing_spreadprec
!!
!>\BRIEF      Spreads the precipitation over the interval chosen based on the interval chosen by the user.
!!
!! DESCRIPTION:	The behaviour of this routine is controlled by the parameter SPREAD_PREC_SEC in the run.def.
!!              The time in second specified by the user will be the one over which the precipitation will last
!!              where the forcing interval has rain or snow.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_spreadprec(time_int, tlen, timebnd_central, time_central, rainf, snowf)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)         ! Time interval to which we will spread precip
    REAL(r_std), INTENT(in)  :: tlen                ! size of time interval in seconds (time step !)
    REAL(r_std), INTENT(in)  :: timebnd_central(:,:)    ! Time interval over which the read data is valid
    REAL(r_std), INTENT(in)  :: time_central(:)     ! Center of the time interval
    REAL(r_std), INTENT(out) :: rainf(:), snowf(:)
    !
    ! LOCAL
    !
    LOGICAL, SAVE :: first_call_spreadprec=.TRUE.
    REAL(r_std), SAVE :: time_to_spread=3600.0
    LOGICAL, SAVE :: spreadprec_cont=.FALSE. 
    !
    INTEGER(i_std) :: imin(1), i, tind(3)
    REAL(r_std) :: ft(3), dt, left, right
    INTEGER(i_std) :: offset, nb_spread
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: tspread
    !
    !
    IF ( .NOT. ALLOCATED(mask) ) THEN
       ALLOCATE(mask(slab_size_max))
       mask(:) = .FALSE.
    ENDIF
    IF ( .NOT. ALLOCATED(tspread) ) THEN 
       ALLOCATE(tspread(nbpoint_proc)) 
       tspread(:) = time_to_spread
    ENDIF
    !
    IF ( first_call_spreadprec ) THEN
       !Config Key   = SPREAD_PREC
       !Config Desc  = On how long we spread the precipitation, value in number of DT_SECHIBA
       !Config Def   = Half of the forcing time step or uniform, depending on dt_force and dt_sechiba
       !Config If    = orchideedriver
       !Config Help  = Spread the precipitation over SPREAD_PREC steps of the splited forcing 
       !Config         time step. This ONLY applied if the forcing time step has been splited.
       !Config         If the value indicated is greater than SPLIT_DT, SPLIT_DT is used for it.
       !Config Units = [-]
       !-
       nb_spread = -1
       CALL getin_p('SPREAD_PREC', nb_spread)
       !
       ! Test if we have read the number of time steps to spread in run.def
       ! If not, then probably the time was given in seconds.
       !
       IF ( nb_spread < 0 ) THEN
          !Config Key   = SPREAD_PREC_SEC
          !Config Desc  = On how long we spread the precipitation, value in seconds
          !Config Def   = Half of the forcing time step
          !Config If    = orchideedriver, only read if SPREAD_PREC is not set
          !Config Help  = Spread the precipitation over n seconds of the forcing time step
          !Config         interval. This ONLY applies when the SPREAD_PREC_SEC is smaller than
          !Config         the forcing time step. Should the user set SPREAD_PREC_SEC=0 we will 
          !Config         assume that the rainfall is uniformely distributed over the forcing interval.
          !Config Units = [seconds]
          
          ! Default value of 'SPREAD_PREC_SEC' if not present in the run.def
          ! if forcing_tstep_ave >= 3 hours: time_to_spread is equal to half of the forcing time step
          ! otherwise, it is equal to the forcing time step 
          IF ( forcing_tstep_ave >= one_day/8.0) THEN
             time_to_spread = forcing_tstep_ave/2.0
          ELSE
             time_to_spread = forcing_tstep_ave
          END IF
             
          CALL getin_p('SPREAD_PREC_SEC', time_to_spread)
       ELSE
          time_to_spread = dt_sechiba_keep * nb_spread
       ENDIF

       !
       ! Add the option to introduce a continuous distribution of precipitation
       !
       !Config Key   = SPREAD_PREC_CONT
       !Config Desc  = Makes the spreading uniform if it still rains at following forcing time step
       !Config Def   = FALSE
       !Config If    = orchideedriver
       !Config Help  = This allows to extend the spreading of rainfall to the entire forcing
       !Config         should it be raining on the following time step. This ensures that if it rains
       !Config         in two consecutive forcing time step at least during the first period rainfall
       !Config         will be continuous. This avoids the spiky nature of rainfall in periods of
       !Config         prolonged rainfall.
       !Config Units = [TRUE/FALSE]
       !
       ! This is the default should 'SPREAD_PREC_CONT' not be present in the run.def
       !
       spreadprec_cont = .FALSE.
       CALL getin_p('SPREAD_PREC_CONT', spreadprec_cont)
       !
       ! Do some verifications on the information read from run.def
       !
       IF ( time_to_spread > forcing_tstep_ave) THEN
          time_to_spread = forcing_tstep_ave
       ELSE IF ( time_to_spread <= 0 ) THEN
          time_to_spread = forcing_tstep_ave
       ENDIF
       !
       first_call_spreadprec = .FALSE.
       !
    ENDIF
    !
    ! First test that we have the right time interval from the forcing to spread the precipitation
    !
    IF ( time_int(1) >= timebnd_central(1,1) .AND. time_int(2) <= timebnd_central(slab_size,2)) THEN
       !
       ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
       !
       mask(1:slab_size) = .TRUE.
       !
       ! To get better precision on the time difference we get a common offset to substract 
       !
       offset = INT(forcingstartdate)
       !
       ! In principle 3 time steps can contribute to the time step closest to the center of the forcing interval
       !
       imin = MINLOC( ABS(time_central(1:slab_size)-(time_int(1)+time_int(2))/2.0), mask(1:slab_size) )
       tind(1) = MAX(imin(1)-1,1)
       tind(2) = imin(1)
       tind(3) = MIN(imin(1)+1,slab_size)
       IF (imin(1)+1 > slab_size) THEN
          WRITE(*,*) "We have a problem here imin(1)+1,slab_size ", imin(1)+1,slab_size
          WRITE(*,*) "Interval : ", time_int(1),time_int(2)
       ENDIF
       !
       tspread(:) = time_to_spread
       !
       DO i=1, nbpoint_proc
          !
          IF ( spreadprec_cont ) THEN
             !
             ! If the next time step has also rainfall, then time_to_spread becomes the length of the forcing time step.
             !
             IF ( rainf_slab(i,tind(3)) > zero .OR. snowf_slab(i,tind(3)) > zero) THEN
                tspread(i) = forcing_tstep_ave
             ELSE
                tspread(i) = time_to_spread
             ENDIF
          ELSE
             ! Default behavious
             tspread(i) = time_to_spread
          ENDIF
          !
          ! Do we need to take some rain from the previous time step ?
          !
          !! Time computation is not better than 1/1000 seconds
          IF ( time_int(1) < timebnd_central(tind(2),1) .AND. preciptime_slab(i,tind(1)) < (tspread(i)-0.001) ) THEN
             dt = ((timebnd_central(tind(2),1)-offset)-(time_int(1)-offset))*one_day
             ft(1) = MIN(tspread(i) - preciptime_slab(i,tind(1)), dt)/tlen
          ELSE
             ft(1) = zero
          ENDIF
          !
          ! Is there still some rain to spread from the current forcing time step ?
          !
          !! Time computation is not better than 1/1000 seconds
          IF (preciptime_slab(i,tind(2)) < (tspread(i)-0.001) ) THEN
             left = MAX(time_int(1), timebnd_central(tind(2),1))
             right = MIN(time_int(2),timebnd_central(tind(2),2))
             dt = ((right-offset)-(left-offset))*one_day
             ft(2) = MIN(tspread(i) - preciptime_slab(i,tind(2)), dt)/tlen
          ELSE
             ft(2) = zero
          ENDIF
          !
          ! Do we need to take some rain from the next time step ?
          !
          !! Time computation is not better than 1/1000 seconds
          IF ( time_int(2) > timebnd_central(tind(2),2) .AND. preciptime_slab(i,tind(3)) < (tspread(i)-0.001) ) THEN
             dt = ((time_int(2)-offset)-(timebnd_central(tind(2),2)-offset))*one_day
             ft(3) = MIN(tspread(i) - preciptime_slab(i,tind(3)), dt)/tlen
          ELSE
             ft(3) = zero
          ENDIF
          !
          ! Do the actual calculation
          !
          rainf(i) = (rainf_slab(i,tind(1)) * forcing_tstep_ave * ft(1) + &
               &  rainf_slab(i,tind(2)) * forcing_tstep_ave * ft(2) + &
               &  rainf_slab(i,tind(3)) * forcing_tstep_ave * ft(3))*tlen/tspread(i)

          snowf(i) = (snowf_slab(i,tind(1)) * forcing_tstep_ave * ft(1) + &
               &  snowf_slab(i,tind(2)) * forcing_tstep_ave * ft(2) + &
               &  snowf_slab(i,tind(3)) * forcing_tstep_ave * ft(3))*tlen/tspread(i)
          !
          ! Update the time over which we have already spread the rainf
          !
          preciptime_slab(i,tind(1)) = preciptime_slab(i,tind(1)) + tlen * ft(1)
          preciptime_slab(i,tind(2)) = preciptime_slab(i,tind(2)) + tlen * ft(2)
          preciptime_slab(i,tind(3)) = preciptime_slab(i,tind(3)) + tlen * ft(3)
          !
       ENDDO
    ELSE
       WRITE(numout,*) "Time interval toward which we will interpolate : ", time_int
       WRITE(numout,*) "Limits of the time slab we have : ", timebnd_central(1,1), timebnd_central(slab_size,2)
       CALL forcing_printdate(time_int(1), "Start of target time interval.")
       CALL forcing_printdate(time_int(2), "End of target time interval.")
       CALL forcing_printdate(timebnd_central(1,1), "Start of time slab we have.")
       CALL forcing_printdate(timebnd_central(slab_size,2), "End of time slab we have.")
       CALL ipslerr (3,'forcing_spreadprec', 'The sitation should not occur Why are we here ?',&
            & "","")
    ENDIF
    

  END SUBROUTINE forcing_spreadprec

!!  =============================================================================================================================
!! SUBROUTINE: forcing_solarint
!!
!>\BRIEF      Interpolates incoming solar radiation to the interval requested.
!!
!! DESCRIPTION:	The interpolation here takes into account the variation of the solar zenith angle
!!              to ensure the diurnal cycle of solar radiation is as well represented as possible.
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_solarint(time_int_in, tlen, timebnd_in, time_cent_in, iim, jjm, lon, lat, model_guess, swdown, solarangle)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)    :: time_int_in(2)             ! Time interval for which we will compute radiation
    REAL(r_std), INTENT(in)    :: tlen                       ! size of time interval in seconds (time step !)
    REAL(r_std), INTENT(in)    :: timebnd_in(:,:)            ! Time interval over which the read data is valid
    REAL(r_std), INTENT(in)    :: time_cent_in(:)            ! Center of the time interval
    INTEGER(i_std), INTENT(in) :: iim, jjm                   ! Size of 2D domain
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm) ! Longitude and latitude
    CHARACTER(LEN=8),INTENT(in) :: model_guess
    REAL(r_std), INTENT(out)   :: swdown(:), solarangle(:)   ! interpolated downward solar radiation and corresponding
    !                                                        ! solar angle.
    !
    ! LOCAL SAVED
    !
    LOGICAL, SAVE        :: first_call_solarint=.TRUE.
    REAL(r_std), SAVE    :: solaryearstart
    INTEGER(i_std), SAVE :: split, split_max
    REAL(r_std), SAVE    :: last_time
    !
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: mean_sinang
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: sinangles
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)     :: time_angles
    ! Dusk-dawn management
    REAL(r_std), SAVE   :: dusk_angle
    !
    ! LOCAL - temporary
    !
    REAL(r_std) :: time_int(2)
    REAL(r_std) :: timebnd(slab_size_max,2)
    REAL(r_std) :: time_cent(slab_size_max)
    INTEGER(i_std) :: year, month, day, hours, minutes
    REAL(r_std) :: sec
    !
    REAL(r_std) :: mean_sol, split_time
    REAL(r_std) :: julian, julian_tmp
    REAL(r_std) :: sinang(iim,jjm)
    INTEGER(i_std) :: is, i, ii, jj, imin(1), tmin(1)
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
    !
    !
    IF ( .NOT. ALLOCATED(mask) ) THEN
       ALLOCATE(mask(slab_size_max))
       mask(:) = .FALSE.
    ENDIF
    
    IF ( first_call_solarint ) THEN
       !
       ! Ensure the offset is on the 1st of Januray of the current years so that we do not
       ! perturbe the solar angle calculations.
       !
       CALL ju2ymds (startdate, year, month, day, sec)
       CALL ymds2ju (year, 1, 1, 0.0, solaryearstart)
       !
       last_time = -9999.0
       !
       ALLOCATE(mean_sinang(iim,jjm))
       mean_sinang(:,:) = 0.0
       !
       split = NINT(forcing_tstep_ave/tlen)
       !
       ! Allow for more space than estimated with the size of the first time step.
       !
       ALLOCATE(time_angles(split*2))
       time_angles(:) = 0.0
       !
       ALLOCATE(sinangles(iim,jjm,split*2))
       sinangles(:,:,:) = 0.0
       !
       dusk_angle=0.01
       !
       first_call_solarint = .FALSE.
       !
       split = 0
       !
    ENDIF
    !
    ! Shift the input dates in order to gain in precision for the time calculations 
    ! solaryearstart: correspond to year, 1, 1,0
    time_int(:) = time_int_in(:)-INT(solaryearstart)
    time_cent(1:slab_size) = time_cent_in(1:slab_size)-INT(solaryearstart)
    timebnd(1:slab_size,1) = timebnd_in(1:slab_size,1)-INT(solaryearstart)
    timebnd(1:slab_size,2) = timebnd_in(1:slab_size,2)-INT(solaryearstart)
    !
    ! Create a mask so that MINLOC does not look outside of the valid interval of time_central
    !
    mask(1:slab_size) = .TRUE.
    !
    ! Locate the time step in the SLAB at hand
    ! find the index in the time slab, which corresponds to the simulation time step
    ! 
    imin = MINLOC( ABS(time_cent(1:slab_size)-(time_int(1)+time_int(2))/2.0), mask(1:slab_size) )
    !
    ! Compute all the angels we will encounter for the current forcing interval
    !
    IF ( last_time .NE. timebnd(imin(1),1) ) THEN
       !
       ! Verify that we have used all the angles of the previous decomposition of the forcing
       ! time step.
       !
       IF ( split .NE. 0 ) THEN
          !
          WRITE(numout,*) "The forcing has a time step of : ", forcing_tstep_ave
          WRITE(numout,*) "The model is configured to run with a time step of : ", tlen
          WRITE(numout,*) "We are left with split = ", split, " starting from ", split_max
          !
          CALL ipslerr (3,'forcing_solarint',"The decomposition of solar downward radiation of the forcing file over the model",&
               &        "has failed. This means the average of the solar radiation over the forcing time step is not conserved.",&
               &        "This can be caused by a time step repeated twice.")
       ENDIF
       !
       ! Compute the number of time steps the model will put in the current interval of forcing data. 
       ! imin: the index corresponding to the  (time_int(1)+time_int(2))/2.0
       ! tmin: the index corresponding to the  (time_int(1)+time_int(2))/2.0, plus 1800 sec
       split = 0
       julian_tmp = (time_int(1)+time_int(2))/2.0
       split_time = julian_tmp+split*tlen/one_day
       tmin = MINLOC( ABS(time_cent(1:slab_size) - split_time), mask(1:slab_size))
       DO WHILE (  tmin(1) .EQ. imin(1) .AND. split_time .LE. timebnd(slab_size,2) )
          split = split + 1
          split_time = julian_tmp+split*tlen/one_day
          tmin = MINLOC( ABS(time_cent(1:slab_size) - split_time), mask(1:slab_size))
       ENDDO
       !
       mean_sinang(:,:) = 0.0
       time_angles(:) = 0.0
       !
       DO is = 1,split
          !
          julian = julian_tmp + (is-1)*tlen/one_day
          !
          ! This call should be better at it allows to compute the difference between the
          ! current date and a reference time to higher precision. But it produces noisy
          ! SWdown fluxes !
!!          CALL solarang (julian, solaryearstart, iim, jjm, lon, lat, sinang)
          ! The old approach.
          ! solarang (iim, jjm)
          IF ( model_guess .EQ. 'unstruct' ) THEN
             CALL solarang_unstruct (julian+INT(solaryearstart), solaryearstart, iim, jjm, lon, lat, sinang)
          ELSE
             CALL solarang (julian+INT(solaryearstart), solaryearstart, iim, jjm, lon, lat, sinang)
          ENDIF
          !
          ! During the dusk,dawn period maintain a minimum angle to take into account the
          ! diffuse radiation which starts before the sun is over the horizon.
          !
          DO jj=1,jjm
             DO ii=1,iim
                IF ( sinang(ii,jj) > zero .AND.  sinang(ii,jj) < dusk_angle ) THEN
                   sinang(ii,jj) = dusk_angle
                ENDIF
                mean_sinang(ii,jj) = mean_sinang(ii,jj)+sinang(ii,jj)
             ENDDO
          ENDDO
          !
          ! Save the solar angle information for later use. That is when we have the target time we will
          ! look in this table the angle we have forseen.
          !
          time_angles(is) = julian
          sinangles(:,:,is) = sinang(:,:)
          !
       ENDDO
       !
       mean_sinang(:,:) = mean_sinang(:,:)/split
       last_time =  timebnd(imin(1),1)
       split_max = split
       !
    ENDIF
    !
    ! For the current timle step get the time of the closest pre-computed solar angle.
    ! 
    julian = (time_int(1)+time_int(2))/2.0
    tmin =  MINLOC(ABS(julian-time_angles(1:split_max)), mask)
    sinang(:,:) = sinangles(:,:,tmin(1))
    ! Remember that we have taken one value of the table for later verification
    split = split - 1
    !
    !orig
    DO i=1, nbpoint_proc
       !
       jj = ((glolindex_proc(i)-1)/iim)+1
       ii = (glolindex_proc(i)-(jj-1)*iim)
       ! swdown_slab: nppoint * timeslab
       IF ( mean_sinang(ii,jj) > zero ) THEN
          swdown(i) = swdown_slab(i,imin(1))*sinang(ii,jj)/mean_sinang(ii,jj)
       ELSE
          swdown(i) = zero
       ENDIF
       !
       ! Why is this ??? Should we not take the angle corresponding to this time step ?? (Jan)
       !
       solarangle(i) = mean_sinang(ii,jj)
       !
    ENDDO    
    !
  END SUBROUTINE forcing_solarint
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_readslab
!!
!>\BRIEF Interface routine to read the data. This routine prepares the memory on each procesor and scatters the read data.
!!
!! DESCRIPTION:
!!
!! \n
!_ ==============================================================================================================================
  !SUBROUTINE forcing_readslab(time_int)
  SUBROUTINE forcing_readslab(time_int, model_guess, is_spinupacc)
    !
    ! This routine serves to interface with forcing_readslab_root and ensure that
    ! the data is distributed correctly on all processors.
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    CHARACTER(LEN=8),INTENT(in) :: model_guess
    LOGICAL, INTENT(in) :: is_spinupacc
    !
    ! Local
    !
    INTEGER(i_std)  :: is, tlen
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: tair_full, qair_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: tairmin_full, tairmax_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: hurs_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: rainf_full, snowf_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: swdown_full, lwdown_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: u_full, v_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: ps_full, ztq_full, zuv_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:)      :: time_precip_old
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: savepreciptime
    INTEGER(i_std) :: s2in1(1)
    LOGICAL, SAVE :: first_call_readslab=.TRUE.
    !
    ! 1.0 Verify that for the slabs the memory is allocated for the variable
    ! as well as its time axis.
    !
    IF ( .NOT. ALLOCATED(tair_slab) ) ALLOCATE(tair_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_tair) ) ALLOCATE(time_tair(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_tair) ) ALLOCATE(timebnd_tair(slab_size,2))
    IF ( forcing_tstep_ave >= one_day/3.0) THEN
       IF ( .NOT. ALLOCATED(tairmax_slab) ) ALLOCATE(tairmax_slab(nbpoint_proc,slab_size))
       IF ( .NOT. ALLOCATED(tairmin_slab) ) ALLOCATE(tairmin_slab(nbpoint_proc,slab_size))
       IF ( .NOT. ALLOCATED(hurs_slab) ) ALLOCATE(hurs_slab(nbpoint_proc,slab_size))
    ENDIF
    !
    IF ( .NOT. ALLOCATED(qair_slab) ) ALLOCATE(qair_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_qair) ) ALLOCATE(time_qair(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_qair) ) ALLOCATE(timebnd_qair(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(rainf_slab) ) ALLOCATE(rainf_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(snowf_slab) ) ALLOCATE(snowf_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_precip) ) ALLOCATE(time_precip(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_precip) ) ALLOCATE(timebnd_precip(slab_size,2))
    IF ( .NOT. ALLOCATED(preciptime_slab) ) ALLOCATE(preciptime_slab(nbpoint_proc,slab_size)) 
    IF ( .NOT. ALLOCATED(time_precip_old) ) ALLOCATE(time_precip_old(slab_size))
    IF ( .NOT. ALLOCATED(savepreciptime) ) ALLOCATE(savepreciptime(nbpoint_proc,slab_size))
    !
    IF ( .NOT. ALLOCATED(swdown_slab) ) ALLOCATE(swdown_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_swdown) ) ALLOCATE(time_swdown(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_swdown) ) ALLOCATE(timebnd_swdown(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(lwdown_slab) ) ALLOCATE(lwdown_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_lwdown) ) ALLOCATE(time_lwdown(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_lwdown) ) ALLOCATE(timebnd_lwdown(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(u_slab) ) ALLOCATE(u_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_u) ) ALLOCATE(time_u(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_u) ) ALLOCATE(timebnd_u(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(v_slab) ) ALLOCATE(v_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_v) ) ALLOCATE(time_v(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_v) ) ALLOCATE(timebnd_v(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(ps_slab) ) ALLOCATE(ps_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(time_ps) ) ALLOCATE(time_ps(slab_size))
    IF ( .NOT. ALLOCATED(timebnd_ps) ) ALLOCATE(timebnd_ps(slab_size,2))
    !
    IF ( .NOT. ALLOCATED(ztq_slab) ) ALLOCATE(ztq_slab(nbpoint_proc,slab_size))
    IF ( .NOT. ALLOCATED(zuv_slab) ) ALLOCATE(zuv_slab(nbpoint_proc,slab_size))    
    !
    ! Save the old time axis for precipitation.
    !
    IF ( first_call_readslab ) THEN
       time_precip_old(:) = zero
       preciptime_slab(:,:) = zero
    ELSE
       time_precip_old(:) = time_precip(:)
    ENDIF
    !
    IF ( is_root_prc) THEN
       !
       ! Allocate on the root processor the memory to receive the data from the file
       !
       ALLOCATE(tair_full(nbpoint_loc,slab_size))
       IF ( forcing_tstep_ave >= one_day/3.0) THEN
          ALLOCATE(tairmax_full(nbpoint_loc,slab_size))
          ALLOCATE(tairmin_full(nbpoint_loc,slab_size))
          ALLOCATE(hurs_full(nbpoint_loc,slab_size))
       ELSE
          ALLOCATE(tairmax_full(1,1))
          ALLOCATE(tairmin_full(1,1))
          ALLOCATE(hurs_full(1,1))
       ENDIF
       ALLOCATE(qair_full(nbpoint_loc,slab_size))
       ALLOCATE(rainf_full(nbpoint_loc,slab_size))
       ALLOCATE(snowf_full(nbpoint_loc,slab_size))
       ALLOCATE(swdown_full(nbpoint_loc,slab_size))
       ALLOCATE(lwdown_full(nbpoint_loc,slab_size))
       ALLOCATE(u_full(nbpoint_loc,slab_size))
       ALLOCATE(v_full(nbpoint_loc,slab_size))
       ALLOCATE(ps_full(nbpoint_loc,slab_size))
       ALLOCATE(ztq_full(nbpoint_loc,slab_size))
       ALLOCATE(zuv_full(nbpoint_loc,slab_size))
       !
       CALL forcing_readslab_root(time_int, first_call_readslab, model_guess, is_spinupacc, &
            &                     tair_full, tairmax_full, tairmin_full, time_tair, timebnd_tair, &
            &                     qair_full, hurs_full, time_qair, timebnd_qair, &
            &                     rainf_full, snowf_full, time_precip, timebnd_precip, &
            &                     swdown_full, time_swdown, timebnd_swdown, &
            &                     lwdown_full, time_lwdown, timebnd_lwdown, &
            &                     u_full, time_u, timebnd_u, &
            &                     v_full, time_v, timebnd_v, &
            &                     ps_full, time_ps, timebnd_ps, &
            &                     ztq_full, zuv_full)
    ELSE
       !
       ALLOCATE(tair_full(1,1))
       IF ( forcing_tstep_ave >= one_day/3.0) THEN
          ALLOCATE(tairmax_full(1,1))
          ALLOCATE(tairmin_full(1,1))
          ALLOCATE(hurs_full(1,1))
       ENDIF
       ALLOCATE(qair_full(1,1))
       ALLOCATE(rainf_full(1,1))
       ALLOCATE(snowf_full(1,1))
       ALLOCATE(swdown_full(1,1))
       ALLOCATE(lwdown_full(1,1))
       ALLOCATE(u_full(1,1))
       ALLOCATE(v_full(1,1))
       ALLOCATE(ps_full(1,1))
       ALLOCATE(ztq_full(1,1))
       ALLOCATE(zuv_full(1,1))
       !
    ENDIF
    !
    ! Broadcast the time information to all procs.
    !
    CALL bcast(slab_size)
    CALL bcast(time_tair)
    CALL bcast(timebnd_tair)
    CALL bcast(time_qair)
    CALL bcast(timebnd_qair)
    CALL bcast(time_precip)
    CALL bcast(timebnd_precip)
    CALL bcast(time_swdown)
    CALL bcast(timebnd_swdown)
    CALL bcast(time_lwdown)
    CALL bcast(timebnd_lwdown)
    CALL bcast(time_u)
    CALL bcast(timebnd_u)
    CALL bcast(time_v)
    CALL bcast(timebnd_v)
    CALL bcast(time_ps)
    CALL bcast(timebnd_ps)
    !
    ! Scatter the slabs of data to all processors
    ! In the classical case, we read and scatter only the land points
    ! but for dynamico, landonly is F. So we do not scatter
    ! If we scatter the data read from dynamico, we read them as 5500*80
    ! then we will scatter all the points (land) to each processor
    
    ! modify the condition here: not only landony but also in the case of dynamico  
    !IF ( landonly  ) THEN
    IF ( landonly .OR. model_guess .EQ. 'unstruct' ) THEN
       CALL scatter(tair_full, tair_slab)
       IF ( forcing_tstep_ave >= one_day/3.0) THEN
          CALL scatter(tairmax_full, tairmax_slab)
          CALL scatter(tairmin_full, tairmin_slab)
          CALL scatter(hurs_full, hurs_slab)
       ENDIF
       CALL scatter(qair_full, qair_slab)
       CALL scatter(rainf_full, rainf_slab)
       CALL scatter(snowf_full, snowf_slab)
       CALL scatter(swdown_full, swdown_slab)
       CALL scatter(lwdown_full, lwdown_slab)
       CALL scatter(u_full, u_slab)
       CALL scatter(v_full, v_slab)
       CALL scatter(ps_full, ps_slab)
       CALL scatter(ztq_full, ztq_slab)
       CALL scatter(zuv_full, zuv_slab)
    ELSE
        tair_slab(:,:) = tair_full(:,:)
       IF ( forcing_tstep_ave >= one_day/3.0) THEN
          tairmax_slab(:,:) = tairmax_full(:,:)
          tairmin_slab(:,:) = tairmin_full(:,:)
          hurs_slab(:,:) = hurs_full(:,:)
       ENDIF
       qair_slab(:,:) = qair_full(:,:)
       rainf_slab(:,:) = rainf_full(:,:)
       snowf_slab(:,:) = snowf_full(:,:)
       swdown_slab(:,:) = swdown_full(:,:)
       lwdown_slab(:,:) = lwdown_full(:,:)
       u_slab(:,:) = u_full(:,:)
       v_slab(:,:) = v_full(:,:)
       ps_slab(:,:) = ps_full(:,:)
       ztq_slab(:,:) = ztq_full(:,:)
       zuv_slab(:,:) = zuv_full(:,:)
    ENDIF
    !
    ! Transfer preciptime_slab from the old slab to the new one.
    ! s2in1: the index of the starting time from new slab in the previous slab
    s2in1 = MINLOC(ABS(time_precip_old-time_precip(1)))
    IF ( size(s2in1) .EQ. 0 ) THEN
       s2in1(1) = 1
    ENDIF
    tlen=slab_size-s2in1(1)+1  
    savepreciptime(:,:) = zero    
    savepreciptime(:,1:tlen) = preciptime_slab(:,s2in1(1):slab_size)    
    preciptime_slab(:,:) = zero
    preciptime_slab(:,1:tlen) = savepreciptime(:,1:tlen)
    
    
    IF (first_call_readslab) THEN
       first_call_readslab = .FALSE.
    ENDIF
    !
    ! Clean-up to free the memory on the root processor.
    !
    IF ( ALLOCATED(tair_full) ) DEALLOCATE(tair_full)
    IF ( ALLOCATED(tairmax_full) ) DEALLOCATE(tairmax_full)
    IF ( ALLOCATED(tairmin_full) ) DEALLOCATE(tairmin_full)
    IF ( ALLOCATED(hurs_full) ) DEALLOCATE(hurs_full)
    IF ( ALLOCATED(qair_full) ) DEALLOCATE(qair_full)
    IF ( ALLOCATED(rainf_full) ) DEALLOCATE(rainf_full)
    IF ( ALLOCATED(snowf_full) ) DEALLOCATE(snowf_full)
    IF ( ALLOCATED(swdown_full) ) DEALLOCATE(swdown_full)
    IF ( ALLOCATED(lwdown_full) ) DEALLOCATE(lwdown_full)
    IF ( ALLOCATED(u_full) ) DEALLOCATE(u_full)
    IF ( ALLOCATED(v_full) ) DEALLOCATE(v_full)
    IF ( ALLOCATED(ps_full) ) DEALLOCATE(ps_full)
    IF ( ALLOCATED(ztq_full) ) DEALLOCATE(ztq_full)
    IF ( ALLOCATED(zuv_full) ) DEALLOCATE(zuv_full)
    !
  END SUBROUTINE forcing_readslab
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_readslab_root
!!
!>\BRIEF Routine which reads a slab of data from the netCDF file and writes it onto the memory.
!!
!! DESCRIPTION:	It is important to read the next slab of data while still keeping an overlap so that
!!              interpolation can continue.
!!              It also attributes a time axis to each variable.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_readslab_root(time_int, first_call_readslab, model_guess, is_spinupacc, &
            &                     tair, tairmax, tairmin, t_tair, tbnd_tair, &
            &                     qair, hurs, t_qair, tbnd_qair, &
            &                     rainf, snowf, t_prec, tbnd_prec, &
            &                     swdown, t_swdown, tbnd_swdown, &
            &                     lwdown, t_lwdown, tbnd_lwdown, &
            &                     u, t_u, tbnd_u, &
            &                     v, t_v, tbnd_v, &
            &                     ps, t_ps, tbnd_ps, &
            &                     ztq, zuv)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    LOGICAL,     INTENT(in)  :: first_call_readslab
    ! added for unstructured grids and spinup tool
    CHARACTER(len=8), INTENT(in) :: model_guess
    LOGICAL, INTENT(in) :: is_spinupacc
    !
    REAL(r_std), INTENT(out) :: tair(:,:), tairmax(:,:), tairmin(:,:)
    REAL(r_std), INTENT(out) :: t_tair(:)
    REAL(r_std), INTENT(out) :: tbnd_tair(:,:)
    !
    REAL(r_std), INTENT(out) :: qair(:,:), hurs(:,:)
    REAL(r_std), INTENT(out) :: t_qair(:)
    REAL(r_std), INTENT(out) :: tbnd_qair(:,:)
    !
    REAL(r_std), INTENT(out) :: rainf(:,:)
    REAL(r_std), INTENT(out) :: snowf(:,:)
    REAL(r_std), INTENT(out) :: t_prec(:)
    REAL(r_std), INTENT(out) :: tbnd_prec(:,:)
    !
    REAL(r_std), INTENT(out) :: swdown(:,:)
    REAL(r_std), INTENT(out) :: t_swdown(:)
    REAL(r_std), INTENT(out) :: tbnd_swdown(:,:)
    !
    REAL(r_std), INTENT(out) :: lwdown(:,:)
    REAL(r_std), INTENT(out) :: t_lwdown(:)
    REAL(r_std), INTENT(out) :: tbnd_lwdown(:,:)
    !
    REAL(r_std), INTENT(out) :: u(:,:)
    REAL(r_std), INTENT(out) :: t_u(:)
    REAL(r_std), INTENT(out) :: tbnd_u(:,:)
    !
    REAL(r_std), INTENT(out) :: v(:,:)
    REAL(r_std), INTENT(out) :: t_v(:)
    REAL(r_std), INTENT(out) :: tbnd_v(:,:)
    !
    REAL(r_std), INTENT(out) :: ps(:,:)
    REAL(r_std), INTENT(out) :: t_ps(:)
    REAL(r_std), INTENT(out) :: tbnd_ps(:,:)
    !
    REAL(r_std), INTENT(out) :: ztq(:,:)
    REAL(r_std), INTENT(out) :: zuv(:,:)
    
    !
    ! Local
    !
    INTEGER(i_std) :: iret, varid
    INTEGER(i_std) :: if, it
    INTEGER(i_std) :: tstart(3), tcount(3)
    INTEGER(i_std) :: imin(1), imax(1), firstif(1)
    INTEGER(i_std) :: nctstart, nctcount, inslabpos
    INTEGER(i_std) :: start_globtime, end_globtime
    INTEGER(i_std) :: timeid_tair, timeid_qair, timeid_precip, timeid_swdown
    INTEGER(i_std) :: timeid_lwdown, timeid_u, timeid_v, timeid_ps, timeid_tmp
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: time_tmp
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: rau
    CHARACTER(LEN=80) :: cellmethod
    
    IF (.NOT. ALLOCATED(time_tmp)) ALLOCATE(time_tmp(slab_size,nbtax))
    IF (.NOT. ALLOCATED(rau)) ALLOCATE(rau(nbpoint_loc,slab_size))
    !
    !
    ! Catch any stupid utilisation of this routine.
    !
    IF ( .NOT. is_root_prc) THEN
       CALL ipslerr (3,'forcing_readslab_root',"Cannot run this routine o other procs than root.",&
            &        "All the information on the forcing files is only lated on the root processor."," ")
    ENDIF
    !
    ! Set some variables to zero
    ! time_int: requested time period for at current simulation time step.
    ! It is in julian dates, with the beginning data in calendar units of forcing file as its beginning date 
    ! time: time array read from forcing file
    ! time(nb_forcing_steps, nbtax*nbtmethods), nbtax*nbtmethods is 12 given the current unstructed grids, 3*4
    ! 
    IF ( first_call_readslab ) THEN
       !
       !
       ! If the first file is only there to provide the last time step of the previous year, we
       ! do not need to read all. We will start 2 forcing time steps before the start of the first
       ! time interval requested.
       ! time: with dimension of total number of data over the three forcing files if there are three.
       ! nbtax: number of time variable or time axis in each forcing file
       ! time(:,1): contains time at the axis 1. In the case of current unstructed grids, we use the time_instant
       ! this is not correct. Because the physical variables uses time_counter (or time_centered)
       ! imin: index in the time variable, which correspond to the value time_int(1)
       imin=MINLOC(ABS(time(:,1)-time_int(1)))
       current_offset = MAX(imin(1)-2,1)
    ELSE      
       !
       ! Compute new offset
       !
       current_offset = position_slab(2)-2    
       !
    ENDIF
    !
    ! Check that the slab size is not too large: this would avoid error if there is only one forcing file.
    !
    IF ( (current_offset-1)+slab_size > nb_forcing_steps) THEN
       slab_size = nb_forcing_steps - (current_offset-1)
    ENDIF
    !
    ! 1.1 Check that the slab we have to read exists
    !
    IF ( slab_size > 0 ) THEN
       !
       start_globtime = current_offset
       end_globtime = (current_offset-1)+slab_size
       inslabpos=1
       WRITE(*,*) ">> Reading from global position ", start_globtime, "up to ", end_globtime
       WRITE(numout,*) ">> Reading from global position ", start_globtime, "up to ", end_globtime
       !
       DO if=MINVAL(time_sourcefile(start_globtime:end_globtime)),MAXVAL(time_sourcefile(start_globtime:end_globtime))
          !
          ! Get position of the part of the global time axis we need to read from this file
          ! means which forcing file 
          firstif = MINLOC(ABS(time_sourcefile-if))
          ! start = distance of start_globtime or nothing + 1 to follow netCDF convention.
          nctstart =  MAX((start_globtime-firstif(1)), 0)+1
          ! count = end index - start index + 1 
          nctcount = MIN((firstif(1)-1)+nbtime_perfile(if),end_globtime)-MAX(firstif(1),start_globtime)+1
          !
          !
          ! Read time over the indexes computed above in order to position the slab correctly
          ! 
          WRITE(*,*) ">> From file ", if," reading from position ", nctstart, "up to ", (nctstart-1)+nctcount
          WRITE(numout,*) ">> From file ", if," reading from position ", nctstart, "up to ", (nctstart-1)+nctcount
          !
          DO it =1,nbtax
             tstart(1) = nctstart
             tcount(1) = nctcount
             ! time_tmp: the seconds since the units of the file
             iret = NF90_GET_VAR(force_id(if), time_id(if,it), time_tmp(inslabpos:inslabpos+nctcount-1,it), tstart, tcount)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(NF90_STRERROR(iret))
                WRITE(*,*) "Working on file ", IF," starting at ",tstart(1)," counting ",tcount(1)
                WRITE(*,*) "The data was to be written in to section ", inslabpos,":",inslabpos+nctcount-1," of time_tmp"
                CALL ipslerr (3,'forcing_readslab_root',"Could not read the time for the current interval."," "," ")
             ENDIF
             ! date0_file: julian date for the beginning (in unit) 
             ! time_tmp: the julian date for the used slab of time
             time_tmp(inslabpos:inslabpos+nctcount-1,it) = date0_file(if,it) + &
                  time_tmp(inslabpos:inslabpos+nctcount-1,it)*convtosec(if)/one_day
             
          ENDDO
          !
          ! 2.0 Find and read variables.
          !
          ! 2.1 Deal with air temperature and humidity as the fist and basic case
          ! hurs: relative humidity
          ! qair: specific humidity
          ! 
          ! orig 202506
          !IF ( ExistVariable(if, "Tair") ) THEN
          ! modified 202506
          IF ( ExistVariable(if, "Tair") .AND. ( forcing_tstep_ave < one_day/3.0)) THEN
             CALL forcing_varforslab(if, "Tair", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, tair, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_tair)
          ELSE IF ( ( forcing_tstep_ave >= one_day/3.0) .AND.  ( ExistVariable(if, "Tmax"))) THEN
             CALL forcing_varforslab(if, "Tmax", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, tairmax, cellmethod)
             CALL forcing_varforslab(if, "Tmin", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, tairmin, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_tair)
             ! orig 202506
             tair(:,:) = tairmax(:,:)
             ! modified 202506
             !IF (ExistVariable(if, "Tair")) THEN
             !   CALL forcing_varforslab(if, "Tair", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, tair, cellmethod)
             !ELSE
             !   tair(:,:) = tairmax(:,:)
             !ENDIF
             
          ELSE
             CALL ipslerr(3, 'forcing_readslab_root','We do not have Tair or Tmax in the file but we need it',&
                     &          'because the timestep of the forcing is sub-diurnal','')
          ENDIF
          
          !
          ! 2.2 Deal with rainfall and snowfall.
          !
          IF ( ExistVariable(if, "Rainf") ) THEN
             CALL forcing_varforslab(if, "Rainf", nctstart, nctcount, inslabpos, model_guess, is_spinupacc,rainf, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_precip)
          
          ELSE IF ( (forcing_tstep_ave >= one_day/3.0)) THEN
             CALL forcing_varforslab(if, "precip", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, rainf, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_precip)
          ENDIF
          !
          IF ( ExistVariable(if, "Snowf")) THEN
             CALL forcing_varforslab(if, "Snowf", nctstart, nctcount, inslabpos,model_guess, is_spinupacc, snowf, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_tmp)
             IF ( timeid_precip .NE. timeid_tmp) THEN
                CALL ipslerr(3, 'forcing_readslab_root','Rainf and Snwof have different time axes.', &
                  &         'Please check the forcing file to ensure both variable have the same cell_method.','')
             ENDIF

          ELSE IF ( (forcing_tstep_ave >= one_day/3.0) )THEN
             snowf(:,:) = rainf(:,:) 
             WHERE (tair(:,:)>=273.15)
                snowf(:,:) = 0.0
             ELSEWHERE (tair(:,:)<273.15)
                rainf(:,:) = 0.0
             END WHERE
          ENDIF
          !
          !
          IF ( ExistVariable(if, "Qair") ) THEN 
             CALL forcing_varforslab(if, "Qair", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, qair, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_qair)
             ! In this case hurs will not be used and can be set to zero.
             hurs(:,:) = zero
          ELSE IF ( ExistVariable(if, "hurs") ) THEN
             CALL forcing_varforslab(if, "hurs", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, hurs, cellmethod)
             CALL forcing_attributetimeaxe(cellmethod, timeid_qair)
             ! In this case qair will not be used and can be set to zero.
             qair(:,:) = zero
             CALL ipslerr(3, 'forcing_readslab_root','hurs is read from forcing file', &
                  &         'Please make sure that the conversion from hurs to specific humidity is correct before continue simulation.','')
          ELSE
             CALL ipslerr(3, 'forcing_readslab_root','We do not have Qair or hurs in the file but we need one of them',&
                     &          'because the timestep of the forcing is sub-diurnal','')
          ENDIF   
          !
          !
          ! 2.3 Deal with downward shorwave and longwave radiation
          !     The SW radiation can have 2 names SWdown_aerosol or SWdown. The first one is
          !     given priority
          !
          CALL forcing_varforslab(if, "SWdown", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, swdown, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_swdown)
          !
          CALL forcing_varforslab(if, "LWdown", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, lwdown, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_lwdown)
          !
          !
          ! 2.4 Deal with wind and pressure
          !
          CALL forcing_varforslab(if, "PSurf", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, ps, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_ps)
          !
          CALL forcing_varforslab(if, "Wind_E", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, u, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_u)
          !
          CALL forcing_varforslab(if, "Wind_N", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, v, cellmethod)
          CALL forcing_attributetimeaxe(cellmethod, timeid_v)
          !
          ! Verify on Tair that we have a credible field.
          ! inslabpos: the indices in spatial domain
          IF (   MINVAL(tair(:,inslabpos:inslabpos+nctcount-1)) < 100.0 .OR. &
               & MAXVAL(tair(:,inslabpos:inslabpos+nctcount-1)) > 500.0 ) THEN
             WRITE(*,*) "ERROR on range of Tair : ", MINVAL(tair(:,inslabpos:inslabpos+nctcount-1)), &
                  &                                  MAXVAL(tair(:,inslabpos:inslabpos+nctcount-1))
             CALL ipslerr(3, 'forcing_readslab_root','The air temperature is not in a credible range.', &
                  &         'Please verify your forcing file.','Are variables for all points to be simulated ?')
          ENDIF
          !
          ! Do the height of the variables T&Q and U&V
          !
          ! First the T&Q level
          !
          IF ( zheight ) THEN
             ztq(:,:) = zlev_fixed
          ELSE IF ( zsigma .OR. zhybrid ) THEN
             DO it=inslabpos,inslabpos+nctcount-1
                rau(:,it) = ps(:,it)/(cte_molr*tair(:,it))
                ztq(:,it) = (ps(:,it) - (zhybrid_a + zhybrid_b*ps(:,it)))/(rau(:,it) * cte_grav)
             ENDDO
          ELSE IF ( zlevels ) THEN
             CALL forcing_varforslab(IF, "Levels", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, ztq, cellmethod)
          ELSE
             CALL ipslerr(3, 'forcing_readslab_root','No case for the vertical levels was specified.', &
                  &         'We cannot determine the height for T and Q.','')
          ENDIF
          !
          ! Now the U&V level
          !
          IF ( zsamelev_uv ) THEN
             zuv(:,:) = ztq(:,:)
          ELSE
             IF ( zheight ) THEN
                zuv(:,:) = zlevuv_fixed
             ELSE IF ( zsigma .OR. zhybrid ) THEN
                DO it=inslabpos,inslabpos+nctcount-1
                   rau(:,it) = ps(:,it)/(cte_molr*tair(:,it))
                   zuv(:,it) = (ps(:,it) - (zhybriduv_a + zhybriduv_b*ps(:,it)))/(rau(:,it) * cte_grav)
                ENDDO
             ELSE IF ( zlevels ) THEN
                CALL forcing_varforslab(IF, "Levels_uv", nctstart, nctcount, inslabpos, model_guess, is_spinupacc, zuv, cellmethod)
             ELSE
                CALL ipslerr(3, 'forcing_readslab_root','No case for the vertical levels was specified.', &
                     &         'We cannot determine the height for U and V.','stop readdim2')
             ENDIF
          ENDIF
          
          inslabpos = inslabpos+nctcount
          
       ENDDO
       !
       ! Use the read time of the slab to place it in the global variables. We consider 
       ! that we can do that on the first axis.
       !
       imin = MINLOC(ABS(time(:,1)-time_tmp(1,1)))
       position_slab(1) = imin(1)
       imax = MINLOC(ABS(time(:,1)-time_tmp(slab_size,1)))
       position_slab(2) = imax(1)
       !
       !
       IF ( position_slab(2)-position_slab(1) .GT. slab_size ) THEN
          DO it =1,nbtax
             WRITE(*,*) "Checking time_tmp on idex : ", it
             WRITE(*,*) "Time_tmp start and end : ",time_tmp(1,it), time_tmp(slab_size,it)
             imin = MINLOC(ABS(time(:,1)-time_tmp(1,it)))
             imax = MINLOC(ABS(time(:,1)-time_tmp(slab_size,it)))
             WRITE(*,*) "Interval read : ", imax(1)-imin(1)+1
          ENDDO
          CALL ipslerr (3,'forcing_readslab_root',"The time slab read does not fit the number of variables read.",&
               &        "Could there be an error in the time axis ?"," ")
       ENDIF
       !
       ! Transfer the global time axis into the time variables approriate for each variable. This way
       ! the time axis for each variable will be centered on the interval of validity. This is an essential assumption 
       ! the time interpolation to be done later.
       !
       WRITE(*,*) "We have found the following axes : ", time_axename(:)
       WRITE(*,*) "For Tair we are using time axis '",TRIM(time_axename(timeid_tair)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_tair)),"."
       t_tair(1:slab_size) = time(position_slab(1):position_slab(2), timeid_tair)
       tbnd_tair(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_tair,:)
       !
       WRITE(*,*) "For Qair we are using time axis '",TRIM(time_axename(timeid_qair)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_qair)),"."
       t_qair(1:slab_size) = time(position_slab(1):position_slab(2), timeid_qair)
       tbnd_qair(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_qair,:)
       !
       WRITE(*,*) "For Rainf and Snowf we are using time axis '",TRIM(time_axename(timeid_precip)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_precip)),"."
       t_prec(1:slab_size) = time(position_slab(1):position_slab(2), timeid_precip)
       tbnd_prec(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_precip,:)
       !
       WRITE(*,*) "For SWdown we are using time axis '",TRIM(time_axename(timeid_swdown)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_swdown)),"."
       t_swdown(1:slab_size) = time(position_slab(1):position_slab(2), timeid_swdown)
       tbnd_swdown(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_swdown,:)
       !
       WRITE(*,*) "For LWdown we are using time axis '",TRIM(time_axename(timeid_lwdown)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_lwdown)),"."
       t_lwdown(1:slab_size) = time(position_slab(1):position_slab(2), timeid_lwdown)
       tbnd_lwdown(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_lwdown,:)
       !
       WRITE(*,*) "For Wind_E we are using time axis '",TRIM(time_axename(timeid_u)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_u)),"."
       t_u(1:slab_size) = time(position_slab(1):position_slab(2), timeid_u)
       tbnd_u(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_u,:)
       !
       WRITE(*,*) "For Wind_N we are using time axis '",TRIM(time_axename(timeid_v)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_v)),"."
       t_v(1:slab_size) = time(position_slab(1):position_slab(2), timeid_v)
       tbnd_v(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_v,:)
       !
       WRITE(*,*) "For PSurf we are using time axis '",TRIM(time_axename(timeid_ps)),&
            &     "' with cell method ",TRIM(time_cellmethod(timeid_ps)),"."
       t_ps(1:slab_size) = time(position_slab(1):position_slab(2), timeid_ps)
       tbnd_ps(1:slab_size,:) = time_bounds(position_slab(1):position_slab(2),timeid_ps,:)
       !
    ELSE
       CALL ipslerr (2,'forcing_readslab_root',"We have reached the end of the slabs we can read.",&
            &          "The calling program needs to manage this situation"," ")
    ENDIF
    !
    ! Have we read to the end of the files ?
    !  
    IF ( current_offset+slab_size >= nb_forcing_steps ) THEN
       end_of_file = .TRUE.
    ELSE
       end_of_file = .FALSE.
    ENDIF
    !
    IF ( ALLOCATED(rau) ) DEALLOCATE(rau)
    IF ( ALLOCATED(time_tmp) ) DEALLOCATE(time_tmp)
    !
  END SUBROUTINE forcing_readslab_root

!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex3d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_reindex3d(nbi, nbj, tin, slab_in, nbout, tout, slab_out, tstart, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj,tin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)
    INTEGER(i_std), INTENT(in) :: tstart
    INTEGER(i_std), INTENT(in) :: reindex(nbout,2)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is, in
    !
    DO is=1,tin
       DO in=1,nbout
          slab_out(in,tstart+(is-1)) = slab_in(reindex(in,1),reindex(in,2),is)
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_reindex3d

!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_reindex2d(nbi, nbj, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout,2)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: in
    !
    DO in=1,nbout
       slab_out(in) = slab_in(reindex(in,1),reindex(in,2))
    ENDDO
    !
  END SUBROUTINE forcing_reindex2d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2dt
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex2dt(nbin, tin, slab_in, nbout, tout, slab_out, tstart, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbin,tin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)
    INTEGER(i_std), INTENT(in) :: tstart
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    ! for compressed crujra halfdeg for example, nbin and nbout are the same, ie., the total land points
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is, in
    !
    DO is=1,tin
       DO in=1,nbout
          slab_out(in,tstart+(is-1)) = slab_in(reindex(in),is)
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcing_reindex2dt

!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2dt_unstruct
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex2dt_unstruct(nbin, tin, slab_in, nbout, tout, slab_out, tstart, lindex, oindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbin,tin)     ! nbin: 16002
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)  ! nbout: 5500
    INTEGER(i_std), INTENT(in) :: tstart
    INTEGER(i_std), INTENT(in) :: lindex(nbout)  ! with the length of 5500, for land points
    INTEGER(i_std), INTENT(in) :: oindex(nbin)
    ! In case of using forcing of dynamico grids:
    ! oindex: with the length of 16002, if landpoint, its values > 0, else it is -1
    !         the value is the indices in referfence to 5500 land points
    ! lindex: with the length of 5500; the values are the indices of land points in reference to 16002 total points
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is, in, ik, ind_total
    !
    DO is=1,tin
       DO in=1,nbin
          ! ik: index in reference to 5500 land points in case of dynamico
          ik = oindex(in)
          ! ind_total: index in refefnce to 16002 points in case of dynamico
          if ( ik > 0) then
             ind_total = lindex(ik)
             slab_out(ik,tstart+(is-1)) = slab_in(ind_total,is)
          endif
       enddo
    enddo
    !
    !
  END SUBROUTINE forcing_reindex2dt_unstruct

  !!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2dt_spinup
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex2dt_spinup(nbin, tin, slab_in, nbout, tout, slab_out, tstart)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, tin, nbout, tout
    REAL(r_std), INTENT(in)    :: slab_in(nbin,tin)     ! nbin: 16002
    REAL(r_std), INTENT(out)   :: slab_out(nbout,tout)  ! nbout: 5500
    INTEGER(i_std), INTENT(in) :: tstart
    !INTEGER(i_std), INTENT(in) :: lindex(nbout)  ! with the length of 5500, for land points
    !INTEGER(i_std), INTENT(in) :: oindex(nbin)
    ! LOCAL
    !
    INTEGER(i_std) :: is, in, ik, ind_total
    ! in the case of spinup tool, the number of input pixels and output pixels are the same, and they are all land pixels
    DO is=1,tin
       DO in=1,nbin
          slab_out(in,tstart+(is-1)) = slab_in(in,is)
       enddo
    enddo
    !
  END SUBROUTINE forcing_reindex2dt_spinup  

!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex1d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex1d(nbin, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbin)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: is
    !
    DO is=1,nbout
       slab_out(is) = slab_in(reindex(is))
    ENDDO
    !
  END SUBROUTINE forcing_reindex1d
!!
!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex2to1
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex2to1(nbi, nbj, slab_in, nbout, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbi, nbj, nbout
    REAL(r_std), INTENT(in)    :: slab_in(nbi,nbj)
    REAL(r_std), INTENT(out)   :: slab_out(nbout)
    INTEGER(i_std), INTENT(in) :: reindex(nbout)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, is
    !
    DO is=1,nbout
       j = INT((reindex(is)-1)/nbi)+1
       i = (reindex(is)-(j-1)*nbi)
       slab_out(is) = slab_in(i,j)
    ENDDO
    !
  END SUBROUTINE forcing_reindex2to1

!!  =============================================================================================================================
!! SUBROUTINE: forcing_reindex1to2
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_reindex1to2(nbin, slab_in, nbi, nbj, slab_out, reindex)
    !
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in) :: nbin, nbi, nbj
    REAL(r_std), INTENT(in)    :: slab_in(nbin)
    REAL(r_std), INTENT(out)   :: slab_out(nbi, nbj)
    INTEGER(i_std), INTENT(in) :: reindex(nbin)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, is
    !
    DO is=1,nbin
       j = INT((reindex(is)-1)/nbi)+1
       i = (reindex(is)-(j-1)*nbi)
       slab_out(i,j) = slab_in(is)
    ENDDO
    !
  END SUBROUTINE forcing_reindex1to2

!!  =============================================================================================================================
!! SUBROUTINE: forcing_close
!!
!>\BRIEF  Close all forcing files
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_close()

    INTEGER(i_std) :: ierr, if

    DO if=1,nb_forcefile
       ierr = NF90_CLOSE(force_id(if))
    ENDDO

  END SUBROUTINE forcing_close

!!  =============================================================================================================================
!! SUBROUTINE: forcing_printdate
!!
!>\BRIEF    Print the date in a human readable format.  
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_printdate(julian_day, message, wunit)
    !
    REAL(r_std), INTENT(in) :: julian_day
    CHARACTER(len=*), INTENT(in) :: message
    INTEGER, INTENT(in), OPTIONAL :: wunit
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    IF (PRESENT(wunit)) THEN
       WRITE(wunit,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2," > ", A60)') &
            &            year, month, day, hours, minutes, seci, message
    ELSE
       WRITE(*,'(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2," > ", A60)') &
            &            year, month, day, hours, minutes, seci, message
    ENDIF
    !
  END SUBROUTINE forcing_printdate

!!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_forgrid
!!
!>\BRIEF     Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_printpoint_forgrid(julian_day, lon_pt, lat_pt, var, message, model_guess)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    REAL(r_std), INTENT(in) :: var(:)
    CHARACTER(len=8), INTENT(in) :: model_guess
    CHARACTER(len=*), INTENT(in) :: message
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    INTEGER(i_std) :: lon_ind, lat_ind, ind
    INTEGER(i_std), DIMENSION(1) :: i,j,k
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       RETURN
    ENDIF
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the local to be analysed
    !
    IF ( model_guess .EQ. 'unstruct') THEN
       i = MINLOC(ABS(lon_loc(1,:)-lon_pt))
       j = MINLOC(ABS(lat_loc(1,:)-lat_pt))
       k = i
       WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F5.1,',', F5.1,'(i=',I6,') Value = ',F12.4,A40)") &
         & hours, minutes, seci, lon_loc(1, i(1)), lat_loc(1,j(1)), k(1), var(k(1)),message
    ELSE
       i = MINLOC(ABS(lon_loc(:,1)-lon_pt))
       j = MINLOC(ABS(lat_loc(1,:)-lat_pt))
       ind = (j(1)-1)*iim_loc+i(1)
       k = MINLOC(ABS(lindex_loc(:)-ind))
       WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F5.1,',', F5.1,'(i=',I6,') Value = ',F12.4,A40)") &
         & hours, minutes, seci, lon_loc(i(1),1), lat_loc(1,j(1)), k(1), var(k(1)), message
    ENDIF
    
    !
    
    
  END SUBROUTINE forcing_printpoint_forgrid
  !!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_forgrid2d
!!
!>\BRIEF     Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_printpoint_forgrid2d(julian_day, lon_pt, lat_pt, var, message)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    REAL(r_std), INTENT(in) :: var(:,:)
    CHARACTER(len=*), INTENT(in) :: message
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec
    INTEGER(i_std) :: lon_ind, lat_ind
    INTEGER(i_std), DIMENSION(1) :: i,j
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       RETURN
    ENDIF
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the local to be analysed
    !
    i = MINLOC(ABS(lon_loc(:,1)-lon_pt))
    j = MINLOC(ABS(lat_loc(1,:)-lat_pt))
    !
    WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F5.1,',', F5.1,'(i=',I6,') Value = ',F12.4,A40)") &
         & hours, minutes, seci, lon_loc(i(1),1), lat_loc(1,j(1)), i(1), j(1), var(i(1),j(1)), message
    
  END SUBROUTINE forcing_printpoint_forgrid2d

!!  =============================================================================================================================
!! SUBROUTINE: forcing_printpoint_gen
!!
!>\BRIEF       Together with the date print some sample values. Useful for checking the forcing. 
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_printpoint_gen(julian_day, lon_pt, lat_pt, nbind, lalo_in, var, message, ktest)
    !
    REAL(r_std), INTENT(in) :: julian_day
    REAL(r_std), INTENT(in) :: lon_pt, lat_pt
    INTEGER(i_std), INTENT(in) :: nbind
    REAL(r_std), INTENT(in) :: lalo_in(:,:)
    REAL(r_std), INTENT(in) :: var(:)
    CHARACTER(len=*), INTENT(in) :: message
    INTEGER(i_std), OPTIONAL, INTENT(out) :: ktest
    !
    !
    !
    INTEGER(i_std) :: year, month, day, hours, minutes, seci
    REAL(r_std) :: sec, mindist
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: dist, refdist
    INTEGER(i_std) :: lon_ind, lat_ind, ind
    INTEGER(i_std) :: i, imin(1)
    REAL(r_std), PARAMETER :: mincos  = 0.0001
    REAL(r_std), PARAMETER :: pi = 3.141592653589793238
    REAL(r_std), PARAMETER :: R_Earth = 6378000.
    !
    ! Check if there is anything to be done
    !
    IF ( MAX(lon_pt, lat_pt) > 360.0 ) THEN
       IF ( PRESENT(ktest) ) ktest = -1
       RETURN
    ENDIF
    !
    ! Allocate memory
    !
    ALLOCATE(dist(nbind))
    ALLOCATE(refdist(nbind))
    !
    ! Convert time first
    !
    CALL ju2ymds (julian_day, year, month, day, sec)
    hours = INT(sec/3600)
    sec = sec - 3600 * hours
    minutes = INT(sec / 60)
    sec = sec - 60 * minutes
    seci = INT(sec)
    !
    ! Get the location to be analysed
    !
    DO i=1,nbind
       dist(i) = acos( sin(lat_pt*pi/180)*sin(lalo_in(i,1)*pi/180) + &
            &    cos(lat_pt*pi/180)*cos(lalo_in(i,1)*pi/180)*&
            &    cos((lalo_in(i,2)-lon_pt)*pi/180) ) * R_Earth
    ENDDO
    !
    ! Look for the next grid point closest to the one with the smalest distance.
    !
    imin = MINLOC(dist)
    DO i=1,nbind
       refdist(i) = acos( sin(lalo_in(imin(1),1)*pi/180)*sin(lalo_in(i,1)*pi/180) + &
            &       cos(lalo_in(imin(1),1)*pi/180)*cos(lalo_in(i,1)*pi/180) * &
            &       cos((lalo_in(i,2)-lalo_in(imin(1),2))*pi/180) ) * R_Earth
    ENDDO
    refdist(imin(1)) =  MAXVAL(refdist)
    mindist = MINVAL(refdist)
    !
    ! Are we closer than the closest points ?
    !
    IF ( PRESENT(ktest) ) ktest = -1
    IF ( dist(imin(1)) <= mindist ) THEN
       !
       WRITE(*,"(I2.2,':',I2.2,':',I2.2,' Loc : ', F6.1,',', F6.1,'(i=',I6,') Value = ',F12.4,A38)") &
            & hours, minutes, seci, lalo_in(imin(1),2), lalo_in(imin(1),1), imin(1), var(imin(1)), message
       !
       IF ( PRESENT(ktest) ) ktest = imin(1)
    ENDIF
    !
  END SUBROUTINE forcing_printpoint_gen

!!  =============================================================================================================================
!! SUBROUTINE: forcing_getglogrid
!!
!>\BRIEF       Routine to read the full grid of the forcing file.
!!
!! DESCRIPTION: The data is stored in the saved variables of the module and serve all other routines. 	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_getglogrid (nbfiles, filename, iim_tmp, jjm_tmp, nbpoint_tmp, closefile, landonly_arg)
    !
    ! This routine reads the global grid information from the forcing file and generates all the 
    ! information needed for this grid.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)    :: nbfiles
    CHARACTER(LEN=*), INTENT(in)  :: filename(:)
    INTEGER(i_std), INTENT(out)   :: iim_tmp, jjm_tmp, nbpoint_tmp
    LOGICAL, INTENT(in)           :: closefile
    LOGICAL, OPTIONAL, INTENT(in) :: landonly_arg
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret, iv, if, lll, nb_glo, ip, ii
    CHARACTER(LEN=20) :: dimname, varname
    CHARACTER(LEN=60) :: lon_units, lat_units, units
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: dimids, londim_ids, latdim_ids
    INTEGER(i_std) :: lon_id, lat_id, land_id, lon_nbdims, lat_nbdims, land_nbdims
    INTEGER(i_std) :: lonvar_id, latvar_id, landvar_id, cell_id, boundslon_id, boundslat_id, area_id
    LOGICAL :: dump_mask, if_unstructured
    INTEGER(i_std) :: ik, i, j, iff, ndimsvar, jj, jp
    ! Read a test variabe
    CHARACTER(len=8) :: testvarname 
    INTEGER(i_std)   :: testvar_id, contfrac_id, iseg
    REAL(r_std) :: testvar_missing, contfrac_missing
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: testvar, contfrac1d, lat1d, lon1d
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: testvar2d, contfrac2d
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: area_full
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: boundslon_full, boundslat_full
    LOGICAL:: is_spinupacc
    REAL(r_std) :: dx, dy, coslat, resolx, resoly
    !
    ! 0.0 Check variables are allocated
    !
    IF ( .NOT. ALLOCATED(force_id)) ALLOCATE(force_id(nbfiles))
    IF ( .NOT. ALLOCATED(id_unlim)) ALLOCATE(id_unlim(nbfiles))
    IF ( .NOT. ALLOCATED(nb_atts)) ALLOCATE(nb_atts(nbfiles))
    IF ( .NOT. ALLOCATED(ndims)) ALLOCATE(ndims(nbfiles))
    IF ( .NOT. ALLOCATED(nvars)) ALLOCATE( nvars(nbfiles))
    !
    ! 0.1 Are we one the root proc ?
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_getglogrid'," This routine can only be called on the root processor.", " ", " ")
    ENDIF
    !
    ! The default behaviour is to provide only land points to the calling program.
    ! But for forcing ocean model there is also the option to pass on land and ocean values.
    ! When the grid is initialized landonly_tmp=.FALSE. has to be set to obtain this behaviour.
    ! orchiee driver scatter only land points on processor. landonly True
    !
    IF ( PRESENT(landonly_arg) ) THEN
       landonly=landonly_arg
    ELSE
       landonly=.TRUE.
    ENDIF
    !
    ! if it is unstructured grids or not
    if_unstructured = .FALSE.
    !
    ! if spinup application is activated or not
    is_spinupacc = .FALSE.
    !
    ! 1.0 Open the netCDF file and get basic dimensions
    !
    DO iff=1,nbfiles
       !
       iret = NF90_OPEN(filename(iff), NF90_NOWRITE, force_id(iff))
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'forcing_getglogrid',"Error opening the forcing file :", filename(iff), " ")
       ENDIF
       !
       iret = NF90_INQUIRE (force_id(iff), nDimensions=ndims(iff), nVariables=nvars(iff), &
            nAttributes=nb_atts(iff), unlimitedDimId=id_unlim(iff))
       !
       !
       ! 2.0 Read the dimensions found in the forcing file. Only deal with the spatial dimension as
       !     time is an unlimited dimension.
       !
       ! initialize nbseg to 4
       nbseg = 4
       !
       DO iv=1,ndims(iff)
          !
          iret = NF90_INQUIRE_DIMENSION (force_id(iff), iv, name=dimname, len=lll)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_getglogrid',"Could not get size of dimensions in file : ", filename(iff), " ")
          ENDIF
          !
          SELECT CASE(lowercase(dimname))
             !
          CASE("west_east")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("south_north")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("longitude")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("latitude")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("lon")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("lat")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("nav_lon")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("nav_lat")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("x")
             CALL forcing_checkdim(iff, filename, iim_glo, lon_id, lll, iv)
          CASE("y")
             CALL forcing_checkdim(iff, filename, jjm_glo, lat_id, lll, iv)
          CASE("land")
             CALL forcing_checkdim(iff, filename, nbland_glo, land_id, lll, iv)
          ! added 20231109 for unstructured grids
          CASE("cell")
             !
             ! If cell dimension exists, this means the unstructured grids used
             ! More computation is realised here in order to find nbland_glo and useful flags
             !
             CALL forcing_checkdim(iff, filename, nb_glo, cell_id, lll, iv)
             iim_glo = 1
             jjm_glo = lll
             !
             ! get the contfrac information
             !
             iret = NF90_INQ_VARID(force_id(1), 'Contfrac', contfrac_id)
             IF (iret /= NF90_NOERR) THEN
                iret = NF90_INQ_VARID(force_id(1), 'contfrac', contfrac_id)
                IF (iret /= NF90_NOERR) THEN
                   CALL ipslerr (3,'forcing_getglogrid',"the variable Contfrac or contfrac not found in the forcing file : ", filename(iff), " ")
                ENDIF
             ENDIF
             
             iret = NF90_GET_ATT(force_id(1), contfrac_id, "_FillValue", contfrac_missing)
             IF (iret /= NF90_NOERR) THEN
                CALL ipslerr (3,'forcing_getglogrid',"the attribute _FillValue is not found for variable Contfrac in the forcing file : ", filename(iff), " ")
             ENDIF
             !
             ! get the number of land point nbland_glo by using contfrac1d values
             !
             IF ( .NOT. ALLOCATED(contfrac1d)) ALLOCATE(contfrac1d(jjm_glo))
             iret = NF90_GET_VAR (force_id(1),contfrac_id,contfrac1d)
             ii = 0
             DO ip=1, jjm_glo
                IF (contfrac1d(ip) > 0. .AND. contfrac1d(ip) < contfrac_missing) THEN
                   ii = ii +1
                ENDIF
             ENDDO
             nbland_glo = ii
             !
             ! If the number of land point nbland_glo is less than total points nb_glo:
             ! then it is unstructured grids with ocean and land included in forcing.
             ! If the number of land point is equal to total points:
             ! then it is the case of spinup accelation and optimization, with only land pixels in the forcing
             ! Redefine the flags landonly and if_unstructured and is_spinup, according to such info
             !
             IF (nbland_glo < nb_glo) THEN
                if_unstructured = .TRUE.
                landonly = .FALSE.
             ELSE IF (nbland_glo .EQ. nb_glo) THEN
                is_spinupacc = .TRUE.
                if_unstructured = .TRUE.
                landonly = .TRUE.
             ENDIF
             !  
          CASE("nvertex")
             !
             ! Currently nvertex is defined in the nc dimension for the unstructured grids only
             ! therefore the relevant variables are allocated below
             ! All the regular grids  
             nbseg = lll
             IF ( .NOT. ALLOCATED(area_full) ) ALLOCATE(area_full(jjm_glo))
             IF ( .NOT. ALLOCATED(boundslon_full) ) ALLOCATE(boundslon_full(nbseg, jjm_glo))
             IF ( .NOT. ALLOCATED(boundslat_full) ) ALLOCATE(boundslat_full(nbseg, jjm_glo))
             IF ( .NOT. ALLOCATED(area_glo)) ALLOCATE(area_glo(iim_glo,jjm_glo))
             IF ( .NOT. ALLOCATED(corners_glo)) ALLOCATE(corners_glo(iim_glo,jjm_glo, nbseg,2))
          END SELECT
          !
       ENDDO
       IF ( iim_glo == 0 .AND. jjm_glo == 0 ) THEN
          CALL ipslerr (3,'forcing_getglogrid',"Did not recognize any dimensions in : ", filename(iff), " ")
       ENDIF
    ENDDO
    !
    ! get right testvarname (often Tair or Tmax) to prepare for the testing of the reading of forcing variable
    !
    iret = NF90_INQ_VARID(force_id(1), 'Tair', iv)
    IF (iret == NF90_NOERR) THEN
       testvarname = 'Tair'
    ELSE
       iret = NF90_INQ_VARID(force_id(1), 'Tmax', iv)
       IF (iret == NF90_NOERR) THEN
          testvarname = 'Tmax'
       ELSE
          CALL ipslerr (3,'forcing_getglogrid',"Could not find variable Tair or Tmax in file."," "," ")
       ENDIF
    ENDIF   
    !
    ! 3.0 Read the spatial coordinate variables found in the first file.
    !
    ALLOCATE(dimids(NF90_MAX_VAR_DIMS), londim_ids(NF90_MAX_VAR_DIMS), latdim_ids(NF90_MAX_VAR_DIMS))
    ! Count the number of time axis we have
    nbtax = 0
    !
    ! initialize variable id in nc file
    contfrac_id = -1
    lonvar_id = -1
    latvar_id = -1
    landvar_id = -1
    testvar_id = -1
    area_id = -1
    boundslat_id = -1
    boundslon_id = -1
    ! nvars(1): 1 represents the first forcing file
    ! nvars(1): represents the max index number of variables in the first forcing file   
    DO iv=1,nvars(1)
       iret = NF90_INQUIRE_VARIABLE(force_id(1), iv, name=varname, ndims=ndimsvar, dimids=dimids)
       iret = NF90_GET_ATT(force_id(1), iv, 'units', units)
       
       ! check if area exist
       IF (  INDEX(lowercase(varname), 'Areas') > 0 .AND. if_unstructured ) THEN 
          iret = NF90_GET_VAR (force_id(1), iv, area_full)    
          ! area_glo: if not land pixel, this variable will have a negative value, not the masked value as areas_full
          DO ip=1,jjm_glo
             area_glo(iim_glo,ip) = area_full(ip)    
          ENDDO
       ENDIF
       
       ! if bounds_lon and lat exist 
       IF ( INDEX(lowercase(varname), "bounds_lon") > 0. .AND. if_unstructured ) THEN
          iret = NF90_GET_VAR (force_id(1), iv, boundslon_full)
          Do iseg=1,nbseg
                DO ip=1,jjm_glo
                ! corners_full should represent the full region or global
                ! corners usually should represent the zoomed region if zoom exists
                ! here is a simplification, by considering the shape of iim jjm, with iim =1 in case of dynamico
                !   corners_glo(iim_glo,jjm_glo, nbseg,2)
                corners_glo(iim_glo,ip,iseg,1) = boundslon_full(iseg, ip)  
             ENDDO
          ENDDO
       ENDIF

       IF ( INDEX(lowercase(varname), "bounds_lat") > 0. .AND. if_unstructured ) THEN
          iret = NF90_GET_VAR (force_id(1), iv, boundslat_full)
          Do iseg=1,nbseg
                DO ip=1,jjm_glo
                ! corners_full should represent the full region or global
                ! corners usually should represent the zoomed region if zoom exists
                ! here is a simplification, by considering the shape of iim jjm, with iim =1
                corners_glo(iim_glo,ip,iseg,2) = boundslat_full(iseg, ip) 
             ENDDO
          ENDDO
       ENDIF
       !
       ! Check that we have the longitude
       !
       IF ( INDEX(lowercase(varname), 'lon') > 0 .AND. INDEX(units, 'east') > 0 ) THEN
          lonvar_id = iv
          lon_units=units
          lon_nbdims = ndimsvar
          londim_ids = dimids
       ENDIF
       !
       ! Check that we have the latitude
       ! INDEX(units, 'east') works, but INDEX(lowercase(units), 'east') gives nothing
       IF ( INDEX(lowercase(varname), 'lat') > 0 .AND. INDEX(units, 'north') > 0) THEN
          latvar_id = iv
          lat_units=units
          lat_nbdims = ndimsvar
          latdim_ids = dimids
       ENDIF
       !
       ! Check that we have the land re-indexing table
       !
       IF ( INDEX(lowercase(varname), 'land') > 0 ) THEN
          landvar_id = iv
          land_nbdims = ndimsvar
          latdim_ids = dimids
       ENDIF
       !
       ! check if we have bounds_lon
       IF ( INDEX(lowercase(varname), "bounds_lon") > 0 ) THEN
          boundslon_id  = iv
       ENDIF
       !
       ! check if we have bounds_lat
       IF ( INDEX(lowercase(varname), "bounds_lat") > 0 ) THEN
          boundslat_id  = iv
       ENDIF
       ! check if we have areas
       !IF ( INDEX(lowercase(varname), "'Areas')") > 0 ) THEN
       IF ( INDEX(lowercase(varname), "Areas") > 0 ) THEN
          area_id  = iv
       ENDIF
       !
       ! Check if we have the contfrac variable
       !
       IF ( INDEX(lowercase(varname), 'contfrac') > 0 ) THEN
          contfrac_id = iv
          iret = NF90_GET_ATT(force_id(1), iv, "_FillValue", contfrac_missing)
          IF (iret /= NF90_NOERR) THEN
             ! No _FillValue found, try to read missing_value instead
             iret = NF90_GET_ATT(force_id(1), iv, "missing_value", contfrac_missing)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(nf90_strerror(iret))
                WRITE(*,*) " >> No _FillValue or missing_value found for contfrac"
                contfrac_missing=1.0E20
             END IF
          ENDIF
       ENDIF
       !
       IF ( INDEX(lowercase(varname), TRIM(lowercase(testvarname))) > 0 .AND. &
            & LEN_TRIM(varname) == LEN_TRIM(testvarname)) THEN
          testvar_id = iv
          iret = NF90_GET_ATT(force_id(1), iv, "_FillValue", testvar_missing)
          IF (iret /= NF90_NOERR) THEN
             ! No _FillValue found, try to read missing_value instead
             iret = NF90_GET_ATT(force_id(1), iv, "missing_value", testvar_missing)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(nf90_strerror(iret))
                WRITE(*,*) " >> No _FillValue or missing_value found for variable=",varname
                testvar_missing=1.0E20
             END IF
          ENDIF
       ENDIF
       !
       ! If we come across time get the calendar and archive it.
       !
       IF ( INDEX(lowercase(units),'seconds since') > 0 .OR. &
            &  INDEX(lowercase(units),'minutes since') > 0 .OR. &
            &  INDEX(lowercase(units),'hours since') > 0 .OR. &
          &  INDEX(lowercase(units),'days since') > 0) THEN 
          !
          ! Get calendar used for the time axis
          !
          iret = NF90_GET_ATT(force_id(1), iv, "calendar", calendar)
          IF (iret == NF90_NOERR) THEN
             WRITE(*,*) ">> Setting the calendar to ",calendar
          ELSE 
             WRITE(*,*) ">> Keeping proleptic Gregorian calendar" 
             calendar="proleptic_gregorian"
          ENDIF
          ! nbtax: number of time axis existing in the forcing file
          ! in the current unstructured grid files, there are 3, time_instant, time_counter and time_centered
          nbtax = nbtax + 1
          !
       ENDIF
    ENDDO
    !
    ! 4.0 Verification that we have found both coordinate variables and the land point index
    !
    IF ( ( lonvar_id < 0 ) .AND. ( INDEX(lon_units, 'east') <= 0 ) ) THEN
       CALL ipslerr (3,'forcing_getglogrid',"Have not found a valid longitude. Units = ", lon_units, " ")
    ENDIF
    !IF ( ( latvar_id < 0 ) .AND. ( INDEX(lowercase(lat_units), 'north') <= 0 ) ) THEN
    IF ( ( latvar_id < 0 ) .AND. ( INDEX(lat_units, 'north') <= 0 ) ) THEN
       CALL ipslerr (3,'forcing_getglogrid',"Have not found a valid latitude. Units = : ", lat_units, " ")
    ENDIF
    IF ( landvar_id < 0 .AND. .NOT. if_unstructured ) THEN
       CALL ipslerr (1,'forcing_getglogrid',"No reindexing table was found. ", &
            &           "This forcing file is not compressed by gathering.", " ")
    ENDIF
    !
    !
    ! 5.0 Allocate the space for the global variables and read them.
    !
    IF ( .NOT. ALLOCATED(lon_glo)) ALLOCATE(lon_glo(iim_glo, jjm_glo))
    IF ( .NOT. ALLOCATED(lat_glo)) ALLOCATE(lat_glo(iim_glo, jjm_glo))
    !
    IF ( lon_nbdims == 2 .AND. lat_nbdims == 2 ) THEN
       iret = NF90_GET_VAR(force_id(1), lonvar_id, lon_glo)
       iret = NF90_GET_VAR(force_id(1), latvar_id, lat_glo)
    ELSE IF ( lon_nbdims == 1 .AND. lat_nbdims == 1 .AND. .NOT. if_unstructured ) THEN
       DO iv=1,jjm_glo
          iret = NF90_GET_VAR(force_id(1), lonvar_id, lon_glo(:,iv))
       ENDDO
       DO iv=1,iim_glo
          iret = NF90_GET_VAR(force_id(1), latvar_id, lat_glo(iv,:))
       ENDDO
    ELSE IF ( lon_nbdims == 1 .AND. lat_nbdims == 1 .AND. if_unstructured ) THEN
       IF ( .NOT. ALLOCATED(lon1d)) ALLOCATE(lon1d(jjm_glo))
       IF ( .NOT. ALLOCATED(lat1d)) ALLOCATE(lat1d(jjm_glo))
       iret = NF90_GET_VAR (force_id(1),lonvar_id,lon1d)
       iret = NF90_GET_VAR (force_id(1),latvar_id,lat1d)
       lon_glo(iim_glo,:) = lon1d(:)
       lat_glo(iim_glo,:) = lat1d(:)
    ELSE
       WRITE(*,*) "Found dimensions for lon and lat :", lon_nbdims, lat_nbdims
       CALL ipslerr (3,'forcing_getglogrid',"Longitude and Latitude variables do not have the right dimensions.", " ", " ")
    ENDIF
    !
    !To find the spatial resolution 
    !
    IF ( if_unstructured .AND. nbseg .EQ. 4 ) THEN
       iret = nf90_get_att(force_id(1), nf90_global, "resolx", resolx)
       IF (iret /= nf90_noerr) THEN
          CALL ipslerr (3,'forcing_getglogrid',"resolx attribute does not exist in the global attribute of forcing file.", " ", " ")
       ENDIF
       iret = nf90_get_att(force_id(1), nf90_global, "resoly", resoly)
       IF (iret /= nf90_noerr) THEN
          CALL ipslerr (3,'forcing_getglogrid',"resoly attribute does not exist in the global attribute of forcing file.", " ", " ")
       ENDIF
    ENDIF
    !
    ! define bound_lat, bound_lon, area and corners in the case of unstructured grids       
    ! if bounds_lon and lat exist 
    IF ( boundslon_id > 0 .AND. if_unstructured ) THEN
       iret = NF90_GET_VAR (force_id(1), boundslon_id, boundslon_full)
       Do iseg=1,nbseg
          DO ip=1,jjm_glo
                ! corners_full should represent the full region or global
                ! corners usually should represent the zoomed region if zoom exists
                ! here is a simplification, by considering the shape of iim jjm, with iim =1 in case of dynamico
                !   corners_glo(iim_glo,jjm_glo, nbseg,2)
             corners_glo(iim_glo,ip,iseg,1) = boundslon_full(iseg, ip)  
          ENDDO
       ENDDO
    ENDIF

    IF ( boundslat_id > 0 .AND. if_unstructured ) THEN
       iret = NF90_GET_VAR (force_id(1), boundslat_id, boundslat_full)
       Do iseg=1,nbseg
          DO ip=1,jjm_glo
             ! corners_full should represent the full region or global
             ! corners usually should represent the zoomed region if zoom exists
             ! here is a simplification, by considering the shape of iim jjm, with iim =1
             corners_glo(iim_glo,ip,iseg,2) = boundslat_full(iseg, ip) 
          ENDDO
       ENDDO
    ENDIF

    IF ( area_id > 0 .AND. if_unstructured ) THEN
       iret = NF90_GET_VAR (force_id(1), area_id, area_full)    
       ! area_glo: if not land pixel, this variable will have a negative value, not the masked value as areas_full
       DO ip=1,jjm_glo
          area_glo(iim_glo,ip) = area_full(ip)    
       ENDDO
    ENDIF
       
    IF ( boundslat_id .LE. 0  .AND. if_unstructured .AND. nbseg .EQ. 4 ) THEN      
       DO j=1,jjm_glo
          corners_glo(1,j,1,1) = lon_glo(1,j)-resolx/2.0
          corners_glo(1,j,2,1) = lon_glo(1,j)+resolx/2.0
          corners_glo(1,j,3,1) = lon_glo(1,j)+resolx/2.0
          corners_glo(1,j,4,1) = lon_glo(1,j)-resolx/2.0
          
          corners_glo(1,j,1,2) = lat_glo(1,j)+resoly/2.0
          corners_glo(1,j,2,2) = lat_glo(1,j)+resoly/2.0
          corners_glo(1,j,3,2) = lat_glo(1,j)-resoly/2.0
          corners_glo(1,j,4,2) = lat_glo(1,j)-resoly/2.0
          
          coslat = MAX( COS(lat_glo(1,j) * pi/180. ), mincos )
          dx = ABS(corners_glo(1,j,2,1) - corners_glo(1,j,1,1)) * pi/180. * R_Earth * coslat
          dy = ABS(corners_glo(1,j,1,2) - corners_glo(1,j,3,2)) * pi/180. * R_Earth
          area_glo(1,j) = dx*dy
          
       ENDDO
    ENDIF
    !
    ! If we have a indexing variable then the data is compressed by gathering, else we have to construct it.
    ! Note: For gloabl regular grids (4 vertex), the compressed flag means that there is a land point dimension
    !       and the forcing variables are indexed along the land points only, so compressed is true in this case
    ! Note: in the case of unstrunctured grids (4 or 6 vertex), there is no land dimension defined in the forcing
    !       But the forcing data are organised in 1D, containing the full domain (land and ocean), or randomly selected points
    !       We set the compressed flag to be True for all the cases of unstructured grids
    ! 
    compressed = .FALSE.
    IF ( landvar_id > 0 ) THEN
       IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbland_glo))
       iret = NF90_GET_VAR(force_id(1), landvar_id, lindex_glo)
       compressed = .TRUE.
    ELSEIF ( landvar_id < 0 .AND. if_unstructured ) THEN
       IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbland_glo))
       compressed = .TRUE.
       jj = 0
       DO jp=1, jjm_glo
          IF (contfrac1d(jp) > 0. .AND. contfrac1d(jp) < contfrac_missing) THEN
             jj = jj + 1
             ! lindex_glo: make the link from land points to their index in the global full domain
             lindex_glo(jj) = jp
          ENDIF
       ENDDO
    ENDIF
    !
    ! Get the land/sea mask and contfrac based on a test variable, if contfract is not available. Else
    ! we use the contfrac variable. The output is to find contfrac_glo, nbpoint_glo, nbland_glo, testvar, lindex_glo
    !
    IF ( .NOT. ALLOCATED(mask_glo)) ALLOCATE(mask_glo(iim_glo, jjm_glo)) 
    IF ( compressed) THEN
       IF ( landonly .OR. if_unstructured ) THEN
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbland_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbland_glo))
          ! total land points (contfrac > 0) : nbland_glo
          ! here testvar is only the land pixels after filtering with contfrac
          ! the output of forcing_contfrac1d is testvar, and contfrac_glo is also defined
          nbpoint_glo = nbland_glo
          CALL forcing_contfrac1d(force_id(1), testvar_id, contfrac_id, testvar, if_unstructured)
       ELSE
          WRITE(*,*) "forcing_tools cannot provide data over ocean points as the"
          WRITE(*,*) "data in the file is compressed by gathering land points."
          WRITE(*,*) "Fatal error"
          CALL ipslerr (3,'forcing_getglogrid',"forcing_tools cannot provide data over ocean points as the", &
               &                               "data in the file is compressed by gathering land points.", " ")
       ENDIF
    ELSE
       IF ( .NOT. ALLOCATED(testvar2d)) ALLOCATE(testvar2d(iim_glo, jjm_glo))
       IF ( .NOT. ALLOCATED(contfrac2d)) ALLOCATE(contfrac2d(iim_glo, jjm_glo))
       CALL forcing_contfrac2d(force_id(1), testvar_id, contfrac_id, testvar2d, contfrac2d, &
            & testvar_missing, contfrac_missing, nbland_glo)
       !
       ! We have found a variable on which we can count the number of land points. We can build
       ! the indexing table and gather the information needed.
       !
       IF ( landonly ) THEN
          nbpoint_glo = nbland_glo
              
          IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbpoint_glo))
          
          IF ( contfrac_id > 0 ) THEN
             ! use the contfrac2d and its missing value to build the land point index
             ! output: 
             ! contfrac_glo: 1D array with contfrac info;
             ! lindex_glo: 1D array with the indices of land point, in reference to the grids in the whole domain
             CALL forcing_buildindex(contfrac2d, contfrac_missing, landonly, lindex_glo, contfrac_glo)
             CALL forcing_reindex(iim_glo, jjm_glo, testvar2d, nbland_glo, testvar, lindex_glo)
          ELSE
             ! use the testvar (often Tair or Tmax) and its missing value to build the land point index
             CALL forcing_buildindex(testvar2d, testvar_missing, landonly, lindex_glo, testvar)
             contfrac_glo(:) = 1.0
             
          ENDIF
       ELSE
          nbpoint_glo = iim_glo*jjm_glo
          IF ( .NOT. ALLOCATED(lindex_glo)) ALLOCATE(lindex_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(contfrac_glo)) ALLOCATE(contfrac_glo(nbpoint_glo))
          IF ( .NOT. ALLOCATED(testvar)) ALLOCATE(testvar(nbpoint_glo))
          IF ( contfrac_id > 0 ) THEN
             CALL forcing_buildindex(contfrac2d, contfrac_missing, landonly, lindex_glo, contfrac_glo)
             CALL forcing_reindex(iim_glo, jjm_glo, testvar2d, nbland_glo, testvar, lindex_glo)
          ELSE
             CALL forcing_buildindex(testvar2d, testvar_missing, landonly, lindex_glo, testvar)
             contfrac_glo(:) = 1.0
          ENDIF
       ENDIF
       !
    ENDIF
    !
    ! We assume that if the forcing file is closed at the end of this subroutine this is a test
    ! or initialisation of the grids. So we will dump the mask in a netCDF file for the user to 
    ! check.
    !
    IF (if_unstructured) THEN
       dump_mask  = .FALSE.
    ELSE
       dump_mask = closefile
    ENDIF
    CALL forcing_checkindex(dump_mask, testvarname, testvar)
    !
    ! 8.0 Close file if needed
    !
    IF ( closefile ) THEN
       CALL forcing_close()
    ENDIF
    !
    ! Prepare variables to be returnned to calling subroutine.
    !
    iim_tmp = iim_glo
    jjm_tmp = jjm_glo
    nbpoint_tmp = nbpoint_glo
    !
    ! Clean up !
    !
    IF ( ALLOCATED(dimids) ) DEALLOCATE(dimids)
    IF ( ALLOCATED(londim_ids) ) DEALLOCATE(londim_ids)
    IF ( ALLOCATED(latdim_ids) ) DEALLOCATE(latdim_ids)
    IF ( ALLOCATED(testvar) ) DEALLOCATE(testvar)
    IF ( ALLOCATED(testvar2d) ) DEALLOCATE(testvar2d)
    IF ( ALLOCATED(contfrac2d) ) DEALLOCATE(contfrac2d)
    IF ( ALLOCATED(contfrac1d) ) DEALLOCATE(contfrac1d)
    !
  END SUBROUTINE forcing_getglogrid

!!  =============================================================================================================================
!! SUBROUTINE: forcing_buildindex
!!
!>\BRIEF      
!!
!! DESCRIPTION:	When the forcing file does not contain compressed variables we need
!!              to build the land index variable from the mask defined by missing variables in 
!!              a test variable.  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_buildindex(var2d, var_missing, landonly, lindex, var)
    !
    ! When the forcing file does not contain compressed variables we need
    ! to build the land index variable from the mask defined by missing variables in 
    ! a test variable.
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in) :: var2d(:,:)
    REAL(r_std), INTENT(in) :: var_missing
    LOGICAL, INTENT(in)     :: landonly
    INTEGER(i_std), INTENT(out) :: lindex(:)
    REAL(r_std), INTENT(out) :: var(:)
    !
    ! Local
    !
    INTEGER(i_std) :: i,j,k
    IF ( MAXVAL(var2d) >= var_missing ) THEN
       ! Case when we have missing values in the vard2d
       k=0
       DO j=1,jjm_glo
          DO i=1,iim_glo
             IF ( landonly ) THEN
                IF ( var2d(i,j) /= var_missing .AND. var2d(i,j) > 0.0 ) THEN
                   k = k + 1
                   lindex(k) = (j-1)*iim_glo+i 
                   var(k) = var2d(i,j)
                ENDIF
             ELSE
                ! When we take all point, no test is performed.
                k = k + 1
                lindex(k) = (j-1)*iim_glo+i 
                var(k) = var2d(i,j)
             ENDIF
          ENDDO
       ENDDO
    ELSE
       ! We suppose that this is land fraction variable
       k=0
       DO j=1,jjm_glo
          DO i=1,iim_glo
             IF ( landonly ) THEN
                IF ( var2d(i,j) > 0.0 ) THEN
                   k = k + 1
                   lindex(k) = (j-1)*iim_glo+i 
                   var(k) = var2d(i,j)
                ENDIF
             ELSE
                ! When we take all point, no test is performed.
                k = k + 1
                lindex(k) = (j-1)*iim_glo+i 
                var(k) = var2d(i,j)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    !
  END SUBROUTINE forcing_buildindex

!!  =============================================================================================================================
!! SUBROUTINE: forcing_contfrac1d
!!
!>\BRIEF      
!!
!! DESCRIPTION:	 This routine build the land/mask if needed and gets the contfrac variable from forcing file.
!!               Here we treat the case where the variables are compressed by gathering. Thus only
!!               land points are available in the file. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_contfrac1d(ifile, testvar_id, contfrac_id, testvar, if_unstructured)
    !
    ! This routine build the land/mask if needed and gets the contfrac variable from forcing file.
    ! Here we treat the case where the variables are compressed by gathering. Thus only
    ! land points are available in the file.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)               :: ifile
    INTEGER(i_std), INTENT(in)               :: testvar_id, contfrac_id
    REAL(r_std), DIMENSION(:), INTENT(inout) :: testvar 
    !
    ! LOCAL
    !
    INTEGER(i_std)                           :: it, iret, ii, ip
    INTEGER(i_std), DIMENSION(3)             :: start, count
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: contfrac2d
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: contfrac1d, testvar1d
    REAL(r_std)                              :: contfrac_missing
    LOGICAL:: if_unstructured
    !
    ! First determine the contfrac variable
    !
    IF ( contfrac_id > 0 ) THEN
       iret = NF90_INQUIRE_VARIABLE(ifile, contfrac_id, ndims=it)
       ! read contfrac and testvar for unstructured grids 
       IF ( if_unstructured ) THEN         
          IF ( .NOT. ALLOCATED(contfrac1d)) ALLOCATE(contfrac1d(jjm_glo))
          IF ( .NOT. ALLOCATED(testvar1d)) ALLOCATE(testvar1d(jjm_glo))
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac1d)
          iret = NF90_GET_ATT(force_id(1), contfrac_id, "_FillValue", contfrac_missing)
          iret = NF90_GET_VAR(ifile, testvar_id, testvar1d)
          ii = 0
          !contfrac_glo: only land pixels
          !testvar: only land pixels
          DO ip=1, jjm_glo
             IF (contfrac1d(ip) > 0. .AND. contfrac1d(ip) < contfrac_missing) THEN
                ii = ii +1
                contfrac_glo(ii) = contfrac1d(ip)
                testvar(ii) = testvar1d(ip)
             ENDIF
          ENDDO
          !
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac1d',"Error reading contfrac for unstructured code ", " ", " ")
          ENDIF

       ELSE IF ( .NOT. if_unstructured .AND. it == 1 ) THEN

          start = (/1,1,0/)
          count = (/nbpoint_glo,1,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac_glo, start, count)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac1d',"Error reading contfrac ", " ", " ")
          ENDIF

       ELSE IF ( it == 2 ) THEN
          IF ( .NOT. ALLOCATED(contfrac2d)) ALLOCATE(contfrac2d(iim_glo,jjm_glo))
          start = (/1,1,0/)
          count = (/iim_glo,jjm_glo,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac2d)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac1d',"Error reading contfrac ", " ", " ")
          ENDIF
          CALL forcing_reindex(iim_glo, jjm_glo, contfrac2d, nbpoint_glo, contfrac_glo, lindex_glo)
          DEALLOCATE(contfrac2d)
       ELSE
          CALL ipslerr (3,'forcing_contfrac1d',"Contfrac has a dimension larger than 2. ", &
               "We do not know how to handle this.", " ")
       ENDIF
    ELSE
       contfrac_glo(:) = 1.0
    ENDIF
    !
    ! Read our test variable 
    !
    iret = NF90_INQUIRE_VARIABLE(ifile, testvar_id, ndims=it)
    ! testvar has dimension of time x space x 
    IF ( it == 2 ) THEN
       ! it=2: time and pixels
       start = (/1,1,0/)
       count = (/nbpoint_glo,1,0/)
    ELSE IF ( it == 3 ) THEN
       ! it =3: time, x, y
       start = (/1,1,1/)
       count = (/nbpoint_glo,1,1/)
       iret = NF90_GET_VAR(ifile, testvar_id, testvar, start, count)
       IF (iret /= NF90_NOERR) THEN
          WRITE(*,*) TRIM(nf90_strerror(iret))
          CALL ipslerr (3,'forcing_contfrac1d',"Error reading testvar.", " ", " ")
       ENDIF
    ELSE
       CALL ipslerr (3,'forcing_contfrac1d',"Testvar has a dimension larger than 3.", &
            "We do not know how to handle this", " ")
    ENDIF
    IF ( ALLOCATED(contfrac1d) ) DEALLOCATE(contfrac1d)
    IF ( ALLOCATED(testvar1d) ) DEALLOCATE(testvar1d)
    IF ( ALLOCATED(contfrac2d) ) DEALLOCATE(contfrac2d)
    !
  END SUBROUTINE forcing_contfrac1d

!!  =============================================================================================================================
!! SUBROUTINE: forcing_contfrac2d
!!
!>\BRIEF      
!!
!! DESCRIPTION: This routine build the land/mask if needed and gets the contfrac variable from forcing file.
!!              Here we treat the case where the variables is 2D. Thus we also need to identify the land points.
!!              For this we can either use the contfrac variable or the test variable.	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_contfrac2d(ifile, testvar_id, contfrac_id, testvar, contfrac, testvar_missing, contfrac_missing, nbland)
    !
    ! This routine build the land/mask if needed and gets the contfrac variable from forcing file.
    ! Here we treat the case where the variables is 2D. Thus we also need to identify the land points.
    ! For this we can either use the contfrac variable or the test variable.
    ! 
    ! ARGUMENTS
    !
    INTEGER(i_std), INTENT(in)                 :: ifile
    INTEGER(i_std), INTENT(in)                 :: testvar_id, contfrac_id
    REAL(r_std), DIMENSION(:,:), INTENT(inout) :: testvar 
    REAL(r_std), DIMENSION(:,:), INTENT(inout) :: contfrac
    REAL(r_std), INTENT(in)                    :: testvar_missing
    REAL(r_std), INTENT(in)                    :: contfrac_missing
    INTEGER(i_std), INTENT(out)                :: nbland
    !
    ! LOCAL
    !
    INTEGER(i_std)                           :: i, j, it, iret
    INTEGER(i_std), DIMENSION(4)             :: start, count
    !
    !
    nbland = 0
    !
    IF ( contfrac_id > 0 ) THEN
       !
       iret = NF90_INQUIRE_VARIABLE(ifile, contfrac_id, ndims=it)
       IF ( it == 2 ) THEN
          start = (/1,1,0,0/)
          count = (/iim_glo,jjm_glo,0,0/)
          iret = NF90_GET_VAR(ifile, contfrac_id, contfrac, start, count)
          IF (iret /= NF90_NOERR) THEN
             WRITE(*,*) TRIM(nf90_strerror(iret))
             CALL ipslerr (3,'forcing_contfrac2d',"Error reading contfrac.", " ", " ")
          ENDIF
       ELSE
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac has a dimension different of 1.", &
               "We do not know how to handle this.", " ")
       ENDIF
       
       IF ( MAXVAL(contfrac) >= contfrac_missing ) THEN
          ! We have missing values in contfrac and we use it to count number of land points
          DO j=1,jjm_glo
             DO i=1,iim_glo
                IF ( contfrac(i,j) /= contfrac_missing .AND. contfrac(i,j) > 0.0 ) THEN
                   nbland = nbland + 1
                ENDIF
             ENDDO
          ENDDO
          
       ELSE
          ! Then ocean is fully contfrc=0 !
          DO j=1,jjm_glo
             DO i=1,iim_glo
                IF ( contfrac(i,j) > 0.0 ) THEN
                   nbland = nbland + 1
                ENDIF
             ENDDO
          ENDDO
          
       ENDIF

       ! If we did not find any land points on the map (i.e. iim_glo > 1 and jjm_glo > 1) then we
       ! look for fractions larger then 0.
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          DO j=1,jjm_glo
             DO i=1,iim_glo
                IF ( contfrac(i,j) > 0.0 ) THEN
                   nbland = nbland + 1
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       
       ! Did we get a result ?
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac was used to count the number of land points.", &
               &          "We still have not found any land points when we looked for contfrac > 0.", " ")
       ENDIF
       !
    ELSE
       ! Just so that we have no un-initialized variable
       contfrac(:,:) = 0.0
    ENDIF
    !
    IF ( testvar_id > 0 ) THEN
       !
       iret = NF90_INQUIRE_VARIABLE(ifile, testvar_id, ndims=it)
       IF ( it == 2 ) THEN
          start = (/1,1,0,0/)
          count = (/iim_glo,jjm_glo,0,0/)
       ELSE IF ( it == 3 ) THEN
          start = (/1,1,1,0/)
          count = (/iim_glo,jjm_glo,1,0/)
       ELSE IF ( it == 4 ) THEN
          start = (/1,1,1,1/)
          count = (/iim_glo,jjm_glo,1,1/)
       ELSE
          CALL ipslerr (3,'forcing_contfrac2d',"testvar has a dimension of 1 or larger than 4.", &
               "We do not know how to handle this.", " ")
       ENDIF
       iret = NF90_GET_VAR(ifile, testvar_id, testvar, start, count)
       IF (iret /= NF90_NOERR) THEN
          WRITE(*,*) TRIM(nf90_strerror(iret))
          CALL ipslerr (3,'forcing_contfrac2d',"Error reading testvar.", " ", " ")
       ENDIF
       !
       ! IF with count frac we did not get the number of land points, let us try it here
       !
       IF ( nbland < 1 ) THEN
          DO j=1,jjm_glo
             DO i=1,iim_glo
                IF ( testvar(i,j) < testvar_missing ) THEN
                   nbland = nbland + 1
                   ! Add infor to contfrac
                   IF ( contfrac_id < 0 ) THEN
                      contfrac(i,j) = 1.0
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
       !
       ! Did we get a result here ?
       !
       IF ( iim_glo > 1 .AND. jjm_glo > 1 .AND. nbland < 1 ) THEN
          CALL ipslerr (3,'forcing_contfrac2d',"Contfrac and testvar were used to count the number", &
               &          "of land points. We have not found any land points.", " ")
       ENDIF
       !
    ENDIF
    !
  END SUBROUTINE forcing_contfrac2d

!!  =============================================================================================================================
!! SUBROUTINE: forcing_checkindex
!!
!>\BRIEF      
!!
!! DESCRIPTION:	 For ORCHIDEE its paralelisation requires that the land points are ordered
!!               in such a way that the longitude runs fastest. That means that we go over the
!!               globle filling one line after the other.
!!               As this might not be the case in a compressed vector of land points, we need to 
!!               put all the points on the 2d grid and then scan them in the right order. 
!!               The reindexing is prepared here.  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_checkindex(dump_mask, testvarname, testvar)
    !
    ! For ORCHIDEE its paralelisation requires that the land points are ordered
    ! in such a way that the longitude runs fastest. That means that we go over the
    ! globle filling one line after the other.
    ! As this might not be the case in a compressed vector of land points, we need to 
    ! put all the points on the 2d grid and then scan them in the right order. 
    ! The reindexing is prepared here.
    !
    LOGICAL          :: dump_mask
    CHARACTER(LEN=*) :: testvarname
    REAL(r_std)      :: testvar(:)
    !
    INTEGER(i_std) :: j, i, ik
    REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: testvar_reind
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:)   :: reindex_check, lindex_check
    !
    !
    !
    ! Get the indices of the land points in the focing file
    !
    !IF ( .NOT. ALLOCATED(reindex_glo)) ALLOCATE(reindex_glo(nbpoint_glo))
    IF ( .NOT. ALLOCATED(origind)) ALLOCATE(origind(iim_glo,jjm_glo))
    IF ( .NOT. ALLOCATED(reindex_check)) ALLOCATE(reindex_check(nbpoint_glo))
    IF ( .NOT. ALLOCATED(lindex_check)) ALLOCATE(lindex_check(nbpoint_glo))
    !
    ! Find the origine of each point in the gathered vector on the xy grid.
    !
    !nbpoint_glo: the number of land pixels only
    !origind:2D array. It makes the connection from 2D global indices (lon, lat) to the 1D land pixels (5500)
    !        the max value of this array is nbpoint_glo
    !lindex_glo: 1D array. It can be read directly from forcing files if land parameter exists.
    !            In the case of unstructured grids, it is created.
    !            It makes the connection from 1D land pixels to the global 1D flatted (lon x lat) indices
    !            the indices in this array are referenced to 16002.
    origind(:,:) = -1
    mask_glo(:,:) = 0
    DO ik=1,nbpoint_glo
       j = INT((lindex_glo(ik)-1)/iim_glo)+1
       i = (lindex_glo(ik)-(j-1)*iim_glo)
       origind(i,j) = ik
       mask_glo(i,j) = 1
    ENDDO
    !
    ! Prepare a reindexing array so that the vector goes in the right order : longitude runs
    ! faster than the latitude. Put then also the right information into lindex_glo.
    !
    ik=0
    DO j=1,jjm_glo
       DO i=1,iim_glo
          IF ( origind(i,j) > zero ) THEN
             ik = ik+1
             reindex_check(ik) = origind(i,j)
             lindex_check(ik) = (j-1)*iim_glo+i
          ENDIF
       ENDDO
    ENDDO
    !
    !
    ! Write the mask and a test variable to a file so that the user can check that all is OK
    !
    IF ( dump_mask) THEN
       !
       ! Scatter the test variable and save it in the file
       !
       WRITE(*,*) MINVAL(testvar), "<< test variable ", TRIM(testvarname), " <<", MAXVAL(testvar) 
       ALLOCATE(testvar_reind(nbpoint_glo))
       !
       CALL forcing_reindex(nbpoint_glo, testvar, nbpoint_glo, testvar_reind, reindex_check)
       !
       
       CALL forcing_writetestvar("forcing_mask_checkglo.nc", iim_glo, jjm_glo, nbpoint_glo, &
            &                    lon_glo(:,1), lat_glo(1,:), lindex_check, mask_glo, &
            &                    testvarname, testvar_reind)
       !
    ENDIF
    !
    ! Clean up !
    !
    IF ( ALLOCATED(testvar_reind) ) DEALLOCATE(testvar_reind)
    IF ( ALLOCATED(reindex_check) ) DEALLOCATE(reindex_check)
    IF ( ALLOCATED(lindex_check) ) DEALLOCATE(lindex_check)
    !
  END SUBROUTINE forcing_checkindex

!!  =============================================================================================================================
!! SUBROUTINE: forcing_writetestvar
!!
!>\BRIEF Write the mask and a test variable to a netCDF file.     
!!
!! DESCRIPTION:	This routine allows to test if the variables read from the forcing files is well read.
!!              Typically the forcing is compressed by gathering and thus it is a safe practice
!!              to verify that the un-compression is done correctly and that all points end-up in the
!!              right place on the global lat/lon grid.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_writetestvar(ncdffile, iim, jjm, nbland, lon, lat, lindex, mask, varname, var)
    !
    ! Write the mask and a test variable to a netCDF file
    !
    ! ARGUMENTS
    !
    CHARACTER(LEN=*), INTENT(in) :: ncdffile
    INTEGER(i_std), INTENT(in)   :: iim, jjm, nbland
    REAL(r_std), INTENT(in)      :: lon(iim), lat(jjm)
    INTEGER(i_std), INTENT(in)   :: lindex(nbland)
    INTEGER(i_std), INTENT(in)   :: mask(iim,jjm)
    CHARACTER(LEN=*), INTENT(in) :: varname
    REAL(r_std), INTENT(in)      :: var(nbland)
    !
    ! Local
    !
    INTEGER(i_std) :: ik, i, j
    INTEGER(i_std) :: iret, nlonid, nlatid, varid, fid, ierr, iland
    INTEGER(i_std) :: testid
    INTEGER(i_std), DIMENSION(2) :: lolaid
    REAL(r_std)    :: test_scat(iim,jjm)
    REAL(r_std)    :: testvar_missing
    !
    !
    test_scat(:,:) = NF90_FILL_REAL
    CALL forcing_reindex(nbland, var, iim, jjm, test_scat, lindex)
    !
    iret = NF90_CREATE(ncdffile, NF90_WRITE, fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'forcing_writetestvar',"Error opening the output file : ", ncdffile, " ")
    ENDIF
    !
    ! Define dimensions
    !
    iret = NF90_DEF_DIM(fid,'lon',iim,lolaid(1))
    iret = NF90_DEF_DIM(fid,'lat',jjm,lolaid(2))
    !
    !
    iret = NF90_DEF_VAR(fid,"lon",NF90_REAL4,lolaid(1),nlonid)
    iret = NF90_PUT_ATT(fid,nlonid,'axis',"X")
    iret = NF90_PUT_ATT(fid,nlonid,'standard_name',"longitude")
    iret = NF90_PUT_ATT(fid,nlonid,'units',"degree_east")
    iret = NF90_PUT_ATT(fid,nlonid,'valid_min',MINVAL(lon_glo))
    iret = NF90_PUT_ATT(fid,nlonid,'valid_max',MAXVAL(lon_glo))
    iret = NF90_PUT_ATT(fid,nlonid,'long_name',"Longitude")
    !
    iret = NF90_DEF_VAR(fid,"lat",NF90_REAL4,lolaid(2),nlatid)
    iret = NF90_PUT_ATT(fid,nlatid,'axis',"Y")
    iret = NF90_PUT_ATT(fid,nlatid,'standard_name',"latitude")
    iret = NF90_PUT_ATT(fid,nlatid,'units',"degree_north")
    iret = NF90_PUT_ATT(fid,nlatid,'valid_min',MINVAL(lat_glo))
    iret = NF90_PUT_ATT(fid,nlatid,'valid_max',MAXVAL(lat_glo))
    iret = NF90_PUT_ATT(fid,nlatid,'long_name',"Latitude")
    !
    iret = NF90_DEF_VAR(fid,"mask",NF90_REAL4,lolaid,varid)
    !
    iret = NF90_DEF_VAR(fid,TRIM(varname),NF90_REAL4,lolaid,testid)
    iret = NF90_PUT_ATT(fid,testid,'_FillValue',NF90_FILL_REAL)
    iret = NF90_PUT_ATT(fid,testid,'missing_value',NF90_FILL_REAL)
    !
    iret = NF90_ENDDEF (fid)
    IF (iret /= NF90_NOERR) THEN
       WRITE(*,*) TRIM(nf90_strerror(iret))
       CALL ipslerr (3,'forcing_writetestvar',"Error ending definitions in file : ", ncdffile, " ")
    ENDIF
    !
    ! Write variables
    !
    iret = NF90_PUT_VAR(fid,nlonid,lon)
    iret = NF90_PUT_VAR(fid,nlatid,lat)
    iret = NF90_PUT_VAR(fid,varid,REAL(mask))
    iret = NF90_PUT_VAR(fid,testid,test_scat)
    !
    ! Close file
    !
    iret = NF90_CLOSE(fid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr (3,'forcing_writetestvar',"Error closing the output file : ", ncdffile, " ")
    ENDIF
    !
  END SUBROUTINE forcing_writetestvar

!!  =============================================================================================================================
!! SUBROUTINE: forcing_zoomgrid 
!!
!>\BRIEF      We zoom into the region requested by the user.
!!
!! DESCRIPTION:	Get the area to be zoomed and the sizes of arrays we will need.
!!              This subroutine fills the *_loc variables.
!!              If requested it will dump a test vraible into a netCDF file.  
!!              *_loc variables does not represent the local processor.
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_zoomgrid (zoom_lon, zoom_lat, filename, model_guess, dumpncdf, iim_jjm_grid, zoom_index_init)
    !
    !
    ! Get the area to be zoomed and the sizes of arrays we will need.
    ! This subroutine fills the *_loc variables.
    ! If requested it will dump a test vraible into a netCDF file.
    ! here _loc means space domain, not the processors on cluster server
    !
    ! ARGUMENTS
    !
    REAL(r_std), DIMENSION(2), INTENT(in) :: zoom_lon, zoom_lat
    CHARACTER(LEN=*), INTENT(in) :: filename , model_guess
    ! in the order of iim, jjm (lon, lat)
    INTEGER, DIMENSION(2), INTENT(in), OPTIONAL :: iim_jjm_grid ! size of zoomed region given by grid file
    !The index of smallest lat lon of the zoomed region, in referene to the full region for the wrf grids
    ! computed in globgrd_getdomsz and pass to this subrouine 
    INTEGER, DIMENSION(2), INTENT(in), OPTIONAL :: zoom_index_init  ! in the order of (lon, lat) 
    LOGICAL, INTENT(in) :: dumpncdf
    REAL(r_std)         :: testvar_missing
    !
    ! LOCAL
    !
    INTEGER(i_std) :: i, j, ic, jc, ik, ig, iseg
    REAL(r_std) :: dx, dy, coslat, resolx, resoly
    REAL(r_std) :: lon_tmp(iim_glo), lat_tmp(jjm_glo)
    REAL(r_std) :: longlo_tmp(iim_glo,jjm_glo)
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_val, lat_val
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: zoom_index
    INTEGER(i_std),DIMENSION(2) ::  tmp_lon_ind1, tmp_lat_ind1, tmp_lon_ind2, tmp_lat_ind2
    INTEGER(i_std) :: iindex_init, jindex_init, iindex_end, jindex_end, i_orig, j_orig
    !
    INTEGER(i_std) :: iret, force_id_zoom, iv, ierr
    INTEGER(i_std), DIMENSION(1) :: imin, jmin
    INTEGER(i_std), DIMENSION(2) :: start, count
    INTEGER(i_std), DIMENSION(3) :: start2d, count2d
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: readvar, zoomedvar
     REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: readvar2d
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: index_glotoloc
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lalo
    CHARACTER(LEN=8) :: testvarname
    !
    ! 0.0 Verify we are on the root processor
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_zoomgrid'," This routine can only be called on the root processor.", " ", " ")
    ENDIF
    !
    ! 0.1 Inform the user
    !
    WRITE(*,*) "Zoom forcing : lon = ", zoom_lon
    WRITE(*,*) "Zoom forcing : lat = ", zoom_lat
    ! Some forcing files have longitudes going from 0 to 360. This code works on the
    ! -180 to 180 longitude grid. So if needed we transform the longitudes of the global grid.
    !
    IF ( MAXVAL(lon_glo) <= 180.0 ) THEN
       longlo_tmp=lon_glo
    ELSE
       DO j=1,jjm_glo
          DO i=1,iim_glo
             IF ( lon_glo(i,j) <= 180.0 ) THEN
                longlo_tmp(i,j) = lon_glo(i,j)
             ELSE
                longlo_tmp(i,j) = lon_glo(i,j)-360
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    !! Define the spatial resolution for the regular grids
    !! resolx and resoly are used to calculate corners later on. 
    IF ( model_guess .EQ. 'regular' ) THEN
       IF (iim_glo==1 .OR. jjm_glo==1) THEN
        ! There is only one pixel or a band of pixels on the forcing file
        resolx=1
        resoly=1
       ELSE
        resolx=ABS(lon_glo(2,1)-lon_glo(1,1))
        resoly=ABS(lat_glo(1,1)-lat_glo(1,2)) 
       ENDIF
    ENDIF
    !
    ! See if we need to zoom
    !
    IF (MINVAL(zoom_lon) > MINVAL(longlo_tmp) .OR. MAXVAL(zoom_lon) < MAXVAL(longlo_tmp) .OR.&
         & MINVAL(zoom_lat) > MINVAL(lat_glo) .OR. MAXVAL(zoom_lat) < MAXVAL(lat_glo) ) THEN
       zoom_forcing = .TRUE.
    ENDIF
    !
    ! Determine the size in x and y of the zoomed region
    ! lon_tmp and lat_tmp: 1 or 0. If 1: in the interested region, if 0: outside the interested region
    ! lon_tmp(iim_glo), lat_tmp(jjm_glo)
    lon_tmp(:) = 0.0
    lat_tmp(:) = 0.0
    IF ( zoom_forcing .AND. model_guess .EQ. 'WRF') THEN
       
       iindex_init = zoom_index_init(1)
       jindex_init = zoom_index_init(2)
       !iindex_end = iindex_init + iim_grid -1
       !jindex_end = jindex_init + jjm_grid -1
       iindex_end = iindex_init + iim_jjm_grid(1) -1
       jindex_end = jindex_init + iim_jjm_grid(2) -1
       
       lat_tmp(jindex_init:jindex_end) = 1.0
       lon_tmp(iindex_init:iindex_end) = 1.0

       iim_loc = iindex_end - iindex_init +1
       jjm_loc = jindex_end - jindex_init +1

       

    ELSE IF ( zoom_forcing  .AND. model_guess .EQ. 'regular') THEN
       
       ! Original code, not modified
       ! Working with the hypothesis it is a regular global grid and bring it back to the -180 to 180 interval
       ! if needed.
       !
       lon_tmp(:) = longlo_tmp(:,1)
       lat_tmp(:) = lat_glo(1,:)
       !
       DO i=1,iim_glo
          IF ( lon_tmp(i) <= MINVAL(zoom_lon) .OR.  lon_tmp(i) >= MAXVAL(zoom_lon) ) THEN
             lon_tmp(i) = 0.0
          ELSE
             lon_tmp(i) = 1.0
          ENDIF
       ENDDO
       DO j=1,jjm_glo
          IF ( lat_tmp(j) <= MINVAL(zoom_lat) .OR. lat_tmp(j) >= MAXVAL(zoom_lat) ) THEN
             lat_tmp(j) = 0.0
          ELSE 
             lat_tmp(j) = 1.0
          ENDIF
       ENDDO
       iim_loc = NINT(SUM(lon_tmp))
       jjm_loc = NINT(SUM(lat_tmp))

    ELSE IF ( .NOT. zoom_forcing  .AND. model_guess .EQ. 'unstruct') THEN
       
       ! works for the full domain of structured or unstructured grids
       ! in case of unstructured grids, lat_tmp has shape 1,
       !
       iindex_init = 1
       jindex_init = 1

       iim_loc = iim_glo
       jjm_loc = jjm_glo
       lon_tmp(:) = 1.0
       lat_tmp(:) = 1.0

    ELSE IF ( .NOT. zoom_forcing  .AND. model_guess .NE. 'unstruct') THEN
       iindex_init = 1
       jindex_init = 1

       iindex_end= size(longlo_tmp(:,1))
       jindex_end= size(lat_glo(1,:))

       iim_loc = iim_glo
       jjm_loc = jjm_glo
       lon_tmp(:) = 1.0
       lat_tmp(:) = 1.0
       
    ENDIF
    !
    ! Determine the number of land points in the zoom
    !
    IF ( .NOT. ALLOCATED(lon_loc) ) ALLOCATE(lon_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(lat_loc) ) ALLOCATE(lat_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(mask_loc) ) ALLOCATE(mask_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(zoom_index) ) ALLOCATE(zoom_index(iim_loc,jjm_loc,2))
    !
    IF ( model_guess .NE. 'unstruct') THEN
       IF ( .NOT. ALLOCATED(lon_val)) ALLOCATE(lon_val(iim_loc))
       IF ( .NOT. ALLOCATED(lat_val)) ALLOCATE(lat_val(jjm_loc))
    ELSE
       IF ( .NOT. ALLOCATED(lon_val)) ALLOCATE(lon_val(jjm_loc))
       IF ( .NOT. ALLOCATED(lat_val)) ALLOCATE(lat_val(jjm_loc))
    ENDIF
    !
    ! Create our new lat/lon system which is in the -180/180 range and South to North and West to East
    ! lon_val, lat_val: get the lat/lon in the interested (zoomed) region
    !
    IF ( model_guess .NE. 'unstruct') THEN
       ic=0
       DO i=1,iim_glo
          IF ( lon_tmp(i) > 0 ) THEN
             ic = ic+1
             lon_val(ic) = longlo_tmp(i,1)
          ENDIF
       ENDDO
       jc=0
       DO j=1,jjm_glo
          IF ( lat_tmp(j) > 0 ) THEN
             jc = jc+1
             lat_val(jc) = lat_glo(1,j)
          ENDIF
       ENDDO
    ELSE
        lon_val(:) = longlo_tmp(1,:)
        lat_val(:) = lat_glo(1,:)
    ENDIF

    IF (ok_sort_lon_lat .AND. model_guess .NE. 'unstruct' ) THEN
       WRITE(*,*) 'Now change order of longitudes and latitudes from the forcing file to increasing order, ok_sort_lon_lat=y' 
       CALL sort(lon_val, iim_loc)
       CALL sort(lat_val, jjm_loc)
    ELSE
       ! unstructured grids have lat/lon in mixed order, not simply organised in increasing or decreasing order.
       ! not to apply the reordering in this case.
       ok_sort_lon_lat=.False.
       WRITE(numout,*) 'No change of the order of longitudes and latitudes from the forcing file, ok_sort_lon_lat=n' 
    END IF

    !
    ! Now find the correspondance between the zoomed & re-ordered grid and the global one.
    !
    IF (  model_guess .EQ. 'unstruct') THEN
       ! 
        ! unstructured grids: only full land domain, no zoomed function (made in 202311)
        ! lat_loc(iim_loc,jjm_loc)), lon_loc(iim_loc,jjm_loc))
        !
        lon_loc(:,:) = longlo_tmp(:,:)
        lat_loc(:,:) = lat_glo(:,:)
        mask_loc(:,:) = mask_glo(:,:)
        !
        ! zoom_index: indices in 3D. Must be consistent with iim and jjm definition.
        zoom_index(1,:,1) = 1 
        zoom_index(1,:,2) = (/(i, i = 1, jjm_glo) /)  
        !
        nbpoint_loc = nbpoint_glo

    ELSE IF ( model_guess .EQ. 'WRF' ) THEN ! .AND.  zoom_forcing) THEN
       ! if WRF grid and zoomed grid, added a few lines here
       mask_loc(:,:) = mask_glo(iindex_init:iindex_end, jindex_init:jindex_end)
       lat_loc(:,:) = lat_glo(iindex_init:iindex_end, jindex_init:jindex_end)
       lon_loc(:,:) = longlo_tmp(iindex_init:iindex_end, jindex_init:jindex_end)
       
       
       nbpoint_loc = 0
       DO j=1,jjm_loc
          DO i=1,iim_loc
             zoom_index(i,j,1) = iindex_init + i - 1
             zoom_index(i,j,2) = jindex_init + j - 1
             ! this is consistent to glogrd_getwrf code
             if ( mask_loc(i,j) .GT. 0.5 ) THEN
                nbpoint_loc = nbpoint_loc +1
             ENDIF
          ENDDO
       ENDDO

    ELSE
       ! if regular grid, or wrf full grids:
       ! the orginal code, not modified
       DO j=1,jjm_loc
          DO i=1,iim_loc
             !
             imin=MINLOC(ABS(longlo_tmp(:,1)-lon_val(i)))
             jmin=MINLOC(ABS(lat_glo(1,:)-lat_val(j)))
             !
             lon_loc(i,j) = longlo_tmp(imin(1),jmin(1))
             lat_loc(i,j) = lat_glo(imin(1),jmin(1))
             mask_loc(i,j) = mask_glo(imin(1),jmin(1))
             !
             zoom_index(i,j,1) = imin(1)
             zoom_index(i,j,2) = jmin(1)
             !
          ENDDO
       ENDDO
       !
       nbpoint_loc = SUM(mask_loc)
       !
    ENDIF
    !
    IF ( .NOT. zoom_forcing .AND. nbpoint_loc .NE. nbpoint_glo) THEN
       WRITE(*,*) "We have not zoomed into the forcing file still we get a number of"
       WRITE(*,*) "land points that is different from what we have in the forcing file."
       STOP "forcing_zoomgrid"
    ENDIF
    !
    IF ( .NOT. ALLOCATED(lindex_loc)) ALLOCATE(lindex_loc(nbpoint_loc))
    IF ( .NOT. ALLOCATED(reindex_loc)) ALLOCATE(reindex_loc(nbpoint_loc))
    IF ( .NOT. ALLOCATED(contfrac_loc)) ALLOCATE(contfrac_loc(nbpoint_loc))
    !
    IF ( .NOT. ALLOCATED(reindex2d_loc)) ALLOCATE(reindex2d_loc(nbpoint_loc,2))
    IF ( .NOT. ALLOCATED(index_glotoloc)) ALLOCATE(index_glotoloc(nbpoint_glo))
    IF ( .NOT. ALLOCATED(lalo)) ALLOCATE(lalo(nbpoint_loc,2))
    !
    ! Do the actual zoom on the grid
    !
    ! Set indices of all points as non existant so that we can fill in as we zoom the 
    ! indices of the points which exist.
    index_glotoloc(:) = -1
    !
    ! Loop only over the zoomed grid
    !
    ! Why does the inner loop need to be ic for the pralalisation ????
    !
    ik = 0
    DO jc=1,jjm_loc
       DO ic=1,iim_loc
          !
          ! Point back from the local to the original global i*j grid
          !
          i = zoom_index(ic,jc,1)
          j = zoom_index(ic,jc,2)
          !
          IF ( origind(i,j) > 0 ) THEN
             ik = ik+1
             ! index of the points in the local grid
             lindex_loc(ik) = (jc-1)*iim_loc+ic
             !
             ! For land points, the index of global grid is saved for the this point on the local grid 
             ! origind represents the index in 1D (compressed in xy, not in x or y)
             ! if there are totally 5500 land points in the case of dynamico, then reindex_loc has the max value of 5500.
             reindex_loc(ik) = origind(i,j)
             !
             ! Keep also the i and j of the global grid for this land point on the local grid
             reindex2d_loc(ik,1) = i
             reindex2d_loc(ik,2) = j
             !
             ! Keep the reverse : on the global grid the location where we put the value of the local grid.
             index_glotoloc(origind(i,j)) = ik
             !
             contfrac_loc(ik) = contfrac_glo(origind(i,j))
             !
             lalo(ik,1) = lat_glo(i,j)
             lalo(ik,2) = longlo_tmp(i,j)
             !
          ENDIF
       ENDDO
    ENDDO
    !
    ! nbland_loc: number of land pixels on global or regional scale
    nbland_loc = 0
    DO ik=1, SIZE(contfrac_loc)
       IF (contfrac_loc(ik) > 0.0) THEN
          nbland_loc = nbland_loc + 1
       ENDIF
    ENDDO
    
    ! ncdfstart, ncdfcount represents the index of pixels in 1D (compressed xy)
    ! only counts for land pixels, ncdfcount is the max number of land points
    ncdfstart = MINVAL(reindex_loc)
    reindex_loc(:) = reindex_loc(:)-ncdfstart+1
    ncdfcount =  MAXVAL(reindex_loc)
    !
    ! Compute the areas and the corners on the grid over which we will run ORCHIDEE.
    ! As this module is dedicated for regular lat/lon forcing we know that we can compute these
    ! variables without further worries.
    !
    IF ( .NOT. ALLOCATED(area_loc)) ALLOCATE(area_loc(iim_loc,jjm_loc))
    IF ( .NOT. ALLOCATED(corners_loc)) ALLOCATE(corners_loc(iim_loc,jjm_loc,nbseg,2))

    IF ( model_guess .EQ. "unstruct" ) THEN
       ! in case of dynamico, Tair(time,16002), not Tair(time, 5500)
       ! this two variables will define the starting and the ending indices in space when reading forcing file
       ncdfstart = 1
       ncdfcount = jjm_loc
    
       ! these _glo are already defined in forcing_getglogrid
       ! need to verify if this is the optimized way of defining these two variables
       area_loc(:,:) = area_glo(:,:)
       DO i=1,jjm_loc
          ! corners usually should represent the zoomed region if zoom exists
          ! here is a simplification, by considering the shape of iim jjm, with iim =1
          corners_loc(iim_loc,i,:,:) = corners_glo(iim_glo,i,:,:)
       ENDDO
      

    ELSE IF ( model_guess .EQ. "WRF" ) THEN
       ! if grid is wrf, with 4 vertex and the same surface area, but with different spatial resolution
       !
       DO j=1,jjm_loc
          DO i=1,iim_loc
             i_orig = iindex_init + i - 1
             j_orig = jindex_init + j - 1
         
             ! a clean proposition, with the first corner on the left and upper, and it turns in clockwise direction.
             CALL grid_tolola(i_orig-0.5, j_orig+0.5, corners_loc(i,j,1,1), corners_loc(i,j,1,2))
             CALL grid_tolola(i_orig+0.5, j_orig+0.5, corners_loc(i,j,2,1), corners_loc(i,j,2,2))
             CALL grid_tolola(i_orig+0.5, j_orig-0.5, corners_loc(i,j,3,1), corners_loc(i,j,3,2))
             CALL grid_tolola(i_orig-0.5, j_orig-0.5, corners_loc(i,j,4,1), corners_loc(i,j,4,2)) 

             ! This is the same as is used in globgrid previously. 
             ! It was put here temporarily to be consistent with previous version of globgrid
             ! but this definition of the first index of grid corner is upper right.
             ! this would produce the inconsistency with dx cy computation in forcing_tools.f90
             !CALL grid_tolola(i_orig+0.5, j_orig+0.5, corners_loc(i,j,1,1), corners_loc(i,j,1,2))
             !CALL grid_tolola(i_orig+0.5, j_orig-0.5, corners_loc(i,j,2,1), corners_loc(i,j,2,2))
             !CALL grid_tolola(i_orig-0.5, j_orig-0.5, corners_loc(i,j,3,1), corners_loc(i,j,3,2))
             !CALL grid_tolola(i_orig-0.5, j_orig+0.5, corners_loc(i,j,4,1), corners_loc(i,j,4,2))

             coslat = MAX( COS(lat_glo(i_orig,j_orig) * pi/180. ), mincos )
             dx = ABS(corners_loc(i,j,2,1) - corners_loc(i,j,1,1)) * pi/180. * R_Earth * coslat
             dy = ABS(corners_loc(i,j,1,2) - corners_loc(i,j,4,2)) * pi/180. * R_Earth
             area_loc(i,j) = dx*dy
          
          ENDDO
       ENDDO

    ELSE
       ! if regular grid with equal degree: original code, not modified
       ! Treat first the longitudes, then treat the latitudes
       !
       DO j=1,jjm_loc
          DO i=1,iim_loc
             corners_loc(i,j,1,1) = lon_loc(i,j)-resolx/2.0
             corners_loc(i,j,2,1) = lon_loc(i,j)+resolx/2.0
             corners_loc(i,j,3,1) = lon_loc(i,j)+resolx/2.0
             corners_loc(i,j,4,1) = lon_loc(i,j)-resolx/2.0
          ENDDO
       ENDDO
       !
       ! Now treat the latitudes
       !
       DO i=1,iim_loc
          DO j=1,jjm_loc
             corners_loc(i,j,1,2) = lat_loc(i,j)+resoly/2.0
             corners_loc(i,j,2,2) = lat_loc(i,j)+resoly/2.0
             corners_loc(i,j,3,2) = lat_loc(i,j)-resoly/2.0
             corners_loc(i,j,4,2) = lat_loc(i,j)-resoly/2.0
          ENDDO
       ENDDO
       !
       ! Compute the areas of the grid (using the simplification that the grid is regular in lon/lat).
       !
       DO i=1,iim_loc
          DO j=1,jjm_loc
             coslat = MAX( COS(lat_loc(i,j) * pi/180. ), mincos )
             dx = ABS(corners_loc(i,j,2,1) - corners_loc(i,j,1,1)) * pi/180. * R_Earth * coslat
             !dy = ABS(corners_loc(i,j,1,2) - corners_loc(i,j,3,2)) * pi/180. * R_Earth
             dy = ABS(corners_loc(i,j,1,2) - corners_loc(i,j,4,2)) * pi/180. * R_Earth
             area_loc(i,j) = dx*dy
          ENDDO
       ENDDO

    ENDIF
    !
    ! If requested we read a variable, zoomin and dump the result into a test file.
    ! dumpncdf is False for unstructured grids by default
    IF ( dumpncdf ) THEN
       iret = NF90_OPEN (filename, NF90_NOWRITE, force_id_zoom)
       IF (iret /= NF90_NOERR) THEN
          CALL ipslerr (3,'forcing_zoomgrid',"Error opening the forcing file :", filename, " ")
       ENDIF
       !
       ALLOCATE(readvar2d(iim_glo,jjm_glo), zoomedvar(nbpoint_loc))
       !
       iret = NF90_INQ_VARID(force_id_zoom, 'Tair', iv)
       IF (iret == NF90_NOERR) THEN
          testvarname = 'Tair'
       ELSE
          iret = NF90_INQ_VARID(force_id_zoom, 'Tmax', iv)
          IF (iret == NF90_NOERR) THEN
             testvarname = 'Tmax'
          ELSE
             CALL ipslerr (3,'forcing_zoomgrid',"Could not find variable Tair or Tmax in file."," "," ")
          ENDIF
       ENDIF
       !
       IF ( compressed .AND. model_guess .NE. 'unstruct') THEN
          !
          start(1) = ncdfstart
          start(2) = 1
          count(1) = ncdfcount
          count(2) = 1
          ALLOCATE(readvar(ncdfcount))
          ! for compressed and regular grids, testvar has the dimension of landpoints (1D)
          ! lat and lon has the dimension of (x, y)
          ! ncdfstart and ncdfcount are 1 and max number of land point 
          ! the following line reads the readvar (with the dimension of land points in case of full domain)
          iret = NF90_GET_VAR(force_id_zoom, iv, readvar, start, count)
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_zoomgrid',"Could not read compressed variable Tair from regular gridded file."," "," ")
          ENDIF
          ! the reindex is simple by calling forcing_reindex1d(nbin, slab_in, nbout, slab_out, reindex)
          CALL forcing_reindex(ncdfcount, readvar, nbpoint_loc, zoomedvar, reindex_loc)
       ELSE
          !
          start2d(1) = 1
          start2d(2) = 1
          start2d(3) = 1
          count2d(1) = iim_glo
          count2d(2) = jjm_glo
          count2d(3) = 1
          !
          iret = NF90_GET_ATT(force_id_zoom, iv, "_FillValue", testvar_missing)
          IF (iret /= NF90_NOERR) THEN
             ! No _FillValue found, try to read missing_value instead
             iret = NF90_GET_ATT(force_id_zoom, iv, "missing_value", testvar_missing)
             IF (iret /= NF90_NOERR) THEN
                WRITE(*,*) TRIM(nf90_strerror(iret))
                testvar_missing=1.0E20
             END IF
          ENDIF
          !
          !
          iret = NF90_GET_VAR(force_id_zoom, iv, readvar2d, start2d, count2d)
          !
          IF (iret /= NF90_NOERR) THEN
             CALL ipslerr (3,'forcing_zoomgrid',"Could not read variable Tair from file."," "," ")
          ENDIF
          CALL forcing_reindex(iim_glo, jjm_glo, readvar2d, nbpoint_loc, zoomedvar, reindex2d_loc)
          !
       ENDIF
       !
       CALL forcing_writetestvar("forcing_mask_zoomloc.nc", iim_loc, jjm_loc, nbpoint_loc, &
            &                    lon_loc(:,1), lat_loc(1,:), lindex_loc, mask_loc, &
            &                    TRIM(testvarname), zoomedvar)
       !
       ierr = NF90_CLOSE(force_id_zoom)
    ENDIF
    !
    ! Clean up
    !
    IF ( ALLOCATED(readvar) ) DEALLOCATE(readvar)
    IF ( ALLOCATED(readvar2d) ) DEALLOCATE(readvar2d)
    IF ( ALLOCATED(zoomedvar) ) DEALLOCATE(zoomedvar)
    IF ( ALLOCATED(index_glotoloc) ) DEALLOCATE(index_glotoloc)
    IF ( ALLOCATED(lalo) ) DEALLOCATE(lalo)
    !
  END SUBROUTINE forcing_zoomgrid

!!  =============================================================================================================================
!! SUBROUTINE: forcing_givegridsize
!!
!>\BRIEF      Routine which exports the size of the grid on which the model will run, i.e. the zoomed grid.
!!
!! DESCRIPTION:	This is needed to transfer the grid information from this module to the glogrid.f90 module.  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_givegridsize (iim, jjm, nblindex)
    !
    ! Provides the size of the land grid to be used to the calling program
    !
    ! Size of the x and y direction of the zoomed area
    INTEGER(i_std), INTENT(out) :: iim, jjm
    ! Number of land points in the zoomed area 
    INTEGER(i_std), INTENT(out) :: nblindex
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_givegridsize'," This routine can only be called on the root processor.", &
            &          "The information requested is only available on root processor.", " ")
    ENDIF
    !
    iim = iim_loc
    jjm = jjm_loc
    nblindex = nbland_loc
    !
  END SUBROUTINE forcing_givegridsize

!!  =============================================================================================================================
!! SUBROUTINE: forcing_
!!
!>\BRIEF      
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_vertical(force_id)
    !
    !- This subroutine explores the forcing file to decide what information is available
    !- on the vertical coordinate.
    !
    INTEGER, INTENT(IN) :: force_id
    !
    INTEGER(i_std) :: iret, ireta, iretb
    !
    INTEGER(i_std) :: sigma_id = -1, sigma_uv_id = -1
    INTEGER(i_std) :: hybsiga_id = -1, hybsiga_uv_id = -1
    INTEGER(i_std) :: hybsigb_id = -1, hybsigb_uv_id = -1
    INTEGER(i_std) :: levels_id = -1, levels_uv_id = -1
    INTEGER(i_std) :: height_id = -1, height_uv_id = -1
    INTEGER(i_std) :: lev_id = -1
    !
    LOGICAL :: var_exists, vara_exists, varb_exists, varuv_exists
    LOGICAL :: foundvar
    LOGICAL :: levlegacy
    !
    !- Set all the defaults
    !
    zfixed=.FALSE.
    zsigma=.FALSE.
    zhybrid=.FALSE.
    zlevels=.FALSE.
    zheight=.FALSE.
    zsamelev_uv = .TRUE.
    levlegacy = .FALSE.
    !
    foundvar = .FALSE.
    !
    !- We have a forcing file to explore so let us see if we find any of the conventions
    !- which allow us to find the height of T,Q,U and V.
    !
    IF ( force_id > 0 ) THEN
       !
       ! Case for sigma levels
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Sigma', sigma_id)
          IF ( (sigma_id >= 0) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zsigma = .TRUE.
             iretb = NF90_INQ_VARID(force_id, 'Sigma_uv', sigma_uv_id)
             IF ( (sigma_uv_id >= 0) .OR. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case for Hybrid levels
       !
       IF ( .NOT. foundvar ) THEN
          var_exists = .FALSE.
          varuv_exists = .FALSE.
          ireta = NF90_INQ_VARID(force_id, 'HybSigA', hybsiga_id)
          IF ( (hybsiga_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             iretb = NF90_INQ_VARID(force_id, 'HybSigB', hybsigb_id)
             IF ( (hybsigb_id >= 0 ) .AND. (iretb == NF90_NOERR) ) THEN
                var_exists=.TRUE.
             ELSE
                CALL ipslerr ( 3, 'forcing_vertical','Missing the B coefficient for', &
                     &         'Hybrid vertical levels for T and Q','stop')
             ENDIF
          ENDIF
          ireta = NF90_INQ_VARID(force_id, 'HybSigA_uv', hybsiga_uv_id)
          IF ( (hybsiga_uv_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             iretb = NF90_INQ_VARID(force_id, 'HybSigB_uv', hybsigb_uv_id)
             IF ( (hybsigb_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) THEN
                varuv_exists=.TRUE.
             ELSE
                CALL ipslerr ( 3, 'forcing_vertical','Missing the B coefficient for', &
                     &         'Hybrid vertical levels for U and V','stop')
             ENDIF
          ENDIF
          IF ( var_exists ) THEN
             foundvar = .TRUE.
             zhybrid = .TRUE.
             IF ( varuv_exists ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case for levels (i.e. a 2d time varying field with the height in meters)
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Levels', levels_id)
          IF ( (levels_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zlevels = .TRUE.
             iretb = NF90_INQ_VARID(force_id, 'Levels_uv', levels_uv_id)
             IF ( (levels_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case where a fixed height is provided in meters
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'Height_Lev1', height_id)
          IF ( (height_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zheight = .TRUE.        
             iretb = NF90_INQ_VARID(force_id, 'Height_Levuv', height_uv_id)
             IF ( (height_uv_id >= 0 ) .AND. (iretb == NF90_NOERR) ) zsamelev_uv = .FALSE.
          ENDIF
       ENDIF
       !
       ! Case where a fixed height is provided in meters in the lev variable
       !
       IF ( .NOT. foundvar ) THEN
          ireta = NF90_INQ_VARID(force_id, 'lev', lev_id)
          IF ( (lev_id >= 0 ) .AND. (ireta == NF90_NOERR) ) THEN
             foundvar = .TRUE.
             zheight = .TRUE.
             levlegacy = .TRUE.
          ENDIF
       ENDIF
       !
    ENDIF
    !
    ! We found forcing variables so we need to extract the values if we are dealing with constant values (i.e. all
    ! except the case zlevels
    !
    IF ( foundvar .AND. .NOT. zlevels ) THEN
       !
       IF ( zheight ) THEN
          !
          ! Constant height
          !
          IF ( levlegacy ) THEN
             iret = NF90_GET_VAR(force_id, lev_id, zlev_fixed)
             IF ( iret /= NF90_NOERR ) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable lev from forcing file in legacy mode', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
          ELSE
             iret = NF90_GET_VAR(force_id, height_id, zlev_fixed)
             IF ( iret /= NF90_NOERR ) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable Height_Lev1 from forcing file', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
             IF ( .NOT. zsamelev_uv ) THEN
                iret = NF90_GET_VAR(force_id, height_uv_id, zlevuv_fixed)
                IF ( iret /= NF90_NOERR ) THEN
                   CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable Height_Levuv from forcing file', &
                        &         'NF90_GET_VAR failed.','stop')
                ENDIF
             ENDIF
          ENDIF
          WRITE(*,*) "forcing_vertical : case ZLEV : Read from forcing file :", zlev_fixed, zlevuv_fixed
          !
       ELSE IF ( zsigma .OR. zhybrid ) THEN
          !
          ! Sigma or hybrid levels
          !
          IF ( zsigma ) THEN
             iret = NF90_GET_VAR(force_id, sigma_id, zhybrid_b)
             zhybrid_a = zero
             IF ( .NOT. zsamelev_uv ) THEN
                iret = NF90_GET_VAR(force_id, sigma_uv_id, zhybriduv_b)
                zhybriduv_a = zero
             ENDIF
          ELSE
             ireta = NF90_GET_VAR(force_id, hybsigb_id, zhybrid_b)
             iretb = NF90_GET_VAR(force_id, hybsiga_id, zhybrid_a)
             IF ( ireta /= NF90_NOERR .OR. iretb /= NF90_NOERR) THEN
                CALL ipslerr ( 3, 'forcing_vertical','Attempted to read variable HybSigA and HybSigB from forcing file', &
                     &         'NF90_GET_VAR failed.','stop')
             ENDIF
             IF ( .NOT. zsamelev_uv ) THEN
                ireta = NF90_GET_VAR(force_id, hybsigb_uv_id, zhybriduv_b)
                iretb = NF90_GET_VAR(force_id, hybsiga_uv_id, zhybriduv_a)
                IF ( ireta /= NF90_NOERR .OR. iretb /= NF90_NOERR) THEN
                   CALL ipslerr ( 3, 'forcing_vertical',&
                        &        'Attempted to read variable HybSigA_uv and HybSigB_uv from forcing file', &
                        &        'NF90_GET_VAR failed.','stop')
                ENDIF
             ENDIF
          ENDIF
          WRITE(*,*) "forcing_vertical : case Pressure coordinates : "
          WRITE(*,*) "Read from forcing file :", zhybrid_b, zhybrid_a, zhybriduv_b, zhybriduv_a
       ELSE
          !
          ! Why are we here ???
          !
          CALL ipslerr ( 3, 'forcing_vertical','What is the option used to describe the height of', &
               &         'the atmospheric forcing ?','Please check your forcing file.')
       ENDIF
    ENDIF
    !
    !- We have no forcing file to explore or we did not find anything. So revert back to the run.def and
    !- read what has been specified by the user.
    !
    IF ( force_id < 0 .OR. .NOT. foundvar ) THEN
       !
       !-
       !Config  Key  = HEIGHT_LEV1
       !Config  Desc = Height at which T and Q are given
       !Config  Def  = 2.0
       !Config  If   = orchideedriver
       !Config  Help = The atmospheric variables (temperature and specific
       !Config         humidity) are measured at a specific level.
       !Config         The height of this level is needed to compute
       !Config         correctly the turbulent transfer coefficients.
       !Config         Look at the description of the forcing
       !Config         DATA for the correct value.
       !Config Units = [m]
       !-
       zlev_fixed = 2.0
       CALL getin('HEIGHT_LEV1', zlev_fixed)
       !-
       !Config  Key  = HEIGHT_LEVW
       !Config  Desc = Height at which the wind is given
       !Config  Def  = 10.0
       !Config  If   = orchideedriver
       !Config  Help = The height at which wind is needed to compute
       !Config         correctly the turbulent transfer coefficients.
       !Config Units = [m]
       !-
       zlevuv_fixed = 10.0
       CALL getin('HEIGHT_LEVW', zlevuv_fixed)

       zheight = .TRUE.

       IF ( ABS(zlevuv_fixed-zlev_fixed) > EPSILON(zlev_fixed)) THEN
          zsamelev_uv = .FALSE.
       ENDIF

       CALL ipslerr ( 2, 'forcing_vertical','The height of the atmospheric forcing variables', &
            &         'was not found in the netCDF file.','Thus the values in run.def were used ... or their defaults.')
    ENDIF

  END SUBROUTINE forcing_vertical

!!  =============================================================================================================================
!! SUBROUTINE: forcing_givegrid
!!
!>\BRIEF      Routine which exports the grid (longitude, latitude, land indices) on which the model will run, i.e. the zoomed grid.
!!
!! DESCRIPTION:	This is needed to transfer the grid information from this module to the glogrid.f90 module.  
!!
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_givegrid (nbseg, lon, lat, mask, area, corners, lindex, contfrac, calendar_tmp)
    !
    ! This subroutine will return to the caller the grid which has been extracted from the
    ! the forcing file. It is assumed that the caller has called forcing_givegridsize before
    ! and knows the dimensions of the fields and thus has done the correct allocations.
    !
    !
    INTEGER(i_std), INTENT(in) :: nbseg
    REAL(r_std), INTENT(out) :: lon(iim_loc,jjm_loc), lat(iim_loc,jjm_loc)
    REAL(r_std), INTENT(out) :: mask(iim_loc,jjm_loc)
    REAL(r_std), INTENT(out) :: area(iim_loc,jjm_loc)
    !REAL(r_std), INTENT(out) :: corners(iim_loc,jjm_loc,4,2)
    REAL(r_std), INTENT(out) :: corners(iim_loc,jjm_loc,nbseg,2)
    INTEGER(i_std), INTENT(out) :: lindex(nbpoint_loc)
    REAL(r_std), INTENT(out) :: contfrac(nbpoint_loc)
    CHARACTER(LEN=20), INTENT(out) :: calendar_tmp
    !
    IF ( .NOT. is_root_prc ) THEN
       CALL ipslerr (3,'forcing_givegrid'," This routine can only be called on the root processor.", &
            &          "The information requested is only available on root processor.", " ")
    ENDIF
    !
    IF (nbpoint_loc .NE. nbland_loc) THEN
       WRITE(numout, *) "forcing_givegrid:: nbpoint_loc=", nbpoint_loc
       WRITE(numout, *) "forcing_givegrid:: nbland_loc=", nbland_loc
       CALL ipslerr(3,'forcing_givegrid','nbpoint_loc and nbland_loc do match', & 
                    'The calculation of land points is not correct','')
    ENDIF
    !
    lon(:,:) = lon_loc(:,:)
    lat(:,:) = lat_loc(:,:)
    !
    mask(:,:) = mask_loc(:,:)
    area(:,:) = area_loc(:,:)
    corners(:,:,:,:) = corners_loc(:,:,:,:)
    !
    !
    lindex(:) = lindex_loc(:)
    contfrac(:) = contfrac_loc(:)
    !
    calendar_tmp = calendar
    !
  END SUBROUTINE forcing_givegrid

!!  =============================================================================================================================
!! SUBROUTINE: forcing_checkdim
!!
!>\BRIEF      
!!
!! DESCRIPTION:	Save the dimension or check that it is equal to the previous value.
!!              Should one of the spatial dimensions be different between 2 files, then we have a big problem.
!!              They probably do not belong to the same set of forcing files.  
!!
!! \n
!_ ==============================================================================================================================

SUBROUTINE forcing_checkdim(ifile, filenames, out_dim, out_id, in_dim, in_id)
  !
  ! Save the dimension or check that it is equal to the previous value.
  ! Should one of the spatial dimensions be different between 2 files, then we have a big problem.
  ! They probably do not belong to the same set of forcing files.
  !
  INTEGER(i_std), INTENT(in) :: ifile
  CHARACTER(LEN=*), INTENT(in) :: filenames(:)
  INTEGER(i_std), INTENT(out) :: out_dim, out_id
  INTEGER(i_std), INTENT(in) :: in_dim, in_id
  !
  IF ( ifile == 1 ) THEN
     out_dim = in_dim
     out_id = in_id
  ELSE
     IF ( out_dim /= in_dim ) THEN
        CALL ipslerr (3,'forcing_ocheckdim', 'The dimension of the file is not the same of the first file opened.', &
             &        'The offending file is : ', filenames(ifile))
     ENDIF
  ENDIF
  !
END SUBROUTINE forcing_checkdim

!!  =============================================================================================================================
!! SUBROUTINE: forcing_time
!!
!>\BRIEF Read the time from each file and create the time axis to be the basis for the simulation.     
!!
!! DESCRIPTION:	This is an important routine which analyses the time axis of the forcing files and
!!              stores the main information in the SAVED variables of this routine. 
!!              As this module manages a list of forcing files we also need to check that the time
!!              axis of all these files is continuous and homogeneous.
!!              The bounds are also build for all the time axes so that we know how to interpret the
!!              various variables. 
!!
!! \n
!_ ==============================================================================================================================

SUBROUTINE forcing_time(nbfiles, filenames)
  !
  ! Read the time from each file and create the time axis to be the basis
  ! for the simulation.
  !
  INTEGER(i_std) :: nbfiles
  CHARACTER(LEN=*) :: filenames(nbfiles)
  !
  INTEGER(i_std) :: iv, it, iff, tcnt, itbase, itbasetmp, ittmp
  INTEGER(i_std) :: tstart, maxtime_infile
  REAL(r_std), ALLOCATABLE, DIMENSION(:) :: timeint, time_read
  REAL(r_std), ALLOCATABLE, DIMENSION(:)   :: time_infiles
  CHARACTER(LEN=20) :: axname, calendar, timevarname
  CHARACTER(LEN=60) :: timestamp, tmpatt, ymd, hms
  INTEGER(i_std) :: tncstart(3), tnccount(3)
  !
  INTEGER(i_std) :: iret, year0, month0, day0, hours0, minutes0, seci
  INTEGER(i_std), DIMENSION(1) :: imax, imin
  REAL(r_std) :: sec0, date_int, date0_tmp, dttmp
  CHARACTER :: strc
  LOGICAL :: check=.FALSE.
  CHARACTER(LEN=20) :: used_calendar, forcing_calendar
  !
  ! Check that we are working on the root proc.
  !
  IF ( .NOT. is_root_prc) THEN
     CALL ipslerr (3,'forcing_time',"Cannot run this routine o other procs than root.",&
          &        "All the information on the forcing files is only lated on the root processor."," ")
  ENDIF
  !
  ! used_calendar is not predefined. the following line will give gregorian which is the default value.
  !
  CALL ioget_calendar(used_calendar)
  !
  ! Size of unlimited dimension added up through the files. If variable not allocated before by another
  ! subroutine, it needs to be done here.
  ! nbtax: number of time axis or time variables in the forcing file.
  ! There are currently 3 in the unstructured codes
  ! date0_file: julian date for the beginning
  IF ( .NOT. ALLOCATED(nbtime_perfile) ) ALLOCATE(nbtime_perfile(nbfiles))
  IF ( .NOT. ALLOCATED(date0_file) ) ALLOCATE(date0_file(nbfiles,nbtax))
  !
  ! Go through all files in the list in order to get the total number of time steps we have
  ! in the nbfiles files to be read.
  ! nb_forcing_steps: total number of forcing time steps over the three full years of forcing
  ! if there are three forcings.
  ! If only one forcign file, then  total number of forcing time steps over one file
  !
  nb_forcing_steps = 0
  maxtime_infile = 0
  DO iff=1,nbfiles
     !
     iret = NF90_INQUIRE_DIMENSION(force_id(iff), id_unlim(iff), name=axname, len=nbtime_perfile(iff))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr (3,'forcing_time',"Could not get size of dimension of unlimited axis"," "," ")
     ENDIF
     nb_forcing_steps =  nb_forcing_steps + nbtime_perfile(iff)
     IF ( nbtime_perfile(iff) > maxtime_infile ) maxtime_infile = nbtime_perfile(iff)
  ENDDO
  !
  ! Allocate the variables needed with the time length just calculated.
  ! These variables are saved in the module
  ! nbtmethods: number of cellmethods possibly for forcing
  ! nb_forcing_steps:total number of forcing time steps over the three full years of forcing
  ALLOCATE(time_infiles(nb_forcing_steps))
  ALLOCATE(time(nb_forcing_steps, nbtax*nbtmethods), time_bounds(nb_forcing_steps,nbtax*nbtmethods,2))
  ALLOCATE(time_axename(nbtax*nbtmethods), time_cellmethod(nbtax*nbtmethods))
  !time_sourcefile: contains the forcing file index for each time step of data
  ALLOCATE(time_sourcefile(nb_forcing_steps))
  ALLOCATE(time_id(nb_forcing_steps, nbtax))
  ! Allocate local variables
  ALLOCATE(time_read(nb_forcing_steps))
  ALLOCATE(timeint(nb_forcing_steps))
  !
  ! Get through all variables to find time_id
  ! The key variables to filled up here are time (the time stamps read in the file) and
  ! time_bounds which give the validity interval for the variables.
  !
  tstart=0
  !
  IF ( check ) WRITE(*,*) "forcing_time : going through ", nbfiles, " files to get the time."
  !
  DO iff=1,nbfiles
     !
     time_id(iff,:)=-1
     !
     ! Go through the variables in the file and find the one which is a time axis.
     !
     tcnt=1
     DO iv=1,nvars(iff)
        iret = NF90_GET_ATT(force_id(iff), iv, "units", tmpatt)
        IF ( INDEX(lowercase(tmpatt),'seconds since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=1.0
        ELSE IF ( INDEX(lowercase(tmpatt),'minutes since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=60.0
        ELSE IF ( INDEX(lowercase(tmpatt),'hours since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=3600.0
        ELSE IF ( INDEX(lowercase(tmpatt),'days since') > 0) THEN
           time_id(iff,tcnt)=iv
           tcnt=tcnt+1
           convtosec(iff)=one_day
        ELSE IF ( INDEX(lowercase(tmpatt),'since') > 0) THEN
           CALL ipslerr ( 2, 'forcing_time', tmpatt, &
            &         'is not the units what we expected (seconds, minutes, hours or days) for time.',&
            &         'It is suggested to use one of the above units in the forcing files.')
        ENDIF
     ENDDO

     IF ( ANY(time_id(iff,:) < 0) ) THEN
        CALL ipslerr (3,'forcing_time',"Incorrect numer of time axes. A time axis is missing ",&
             &        "in file :", filenames(iff))
     ENDIF
     !
     IF ( check ) WRITE(*,*) "forcing_time : Looking at time axis for file ", force_id(iff)
     !
     ! Looping through the time axes (or time variables) and read them.
     !
     DO tcnt=1,nbtax
        ! timevarname: the name of time variable in the forcing file
        ! In the current unstructured grids, there are three time variables, time_instant, time_counter and time_centered
        iret = NF90_INQUIRE_VARIABLE(force_id(iff), time_id(iff,tcnt), name=timevarname)
        IF ( check ) WRITE(*,*) "forcing_time : in ", iff, " found variable ", timevarname
        WRITE(numout,*) "forcing_time : in ", iff, " found variable ", timevarname
        !
        ! Get units of time axis
        !
        iret = NF90_GET_ATT(force_id(iff), time_id(iff,tcnt), "units", timestamp) 
        IF ( check ) WRITE(*,*) "forcing_time : has time stamp ", timestamp
        !
        ! Transform the start date of the netCDF file into a julian date for the model
        !
        timestamp = TRIM(timestamp(INDEX(timestamp,'since')+6:LEN_TRIM(timestamp)))
        ymd=TRIM(timestamp(1:INDEX(timestamp,' ')))
        hms=TRIM(timestamp(INDEX(timestamp,' ')+1:LEN_TRIM(timestamp)))
        !
        ! First extral the information from the date string
        !
        READ(ymd(1:INDEX(ymd,'-')-1),'(I4)') year0
        ymd=TRIM(ymd(INDEX(ymd,'-')+1:LEN_TRIM(ymd)))
        READ(ymd(1:INDEX(ymd,'-')-1),'(I2)') month0
        ymd=TRIM(ymd(INDEX(ymd,'-')+1:LEN_TRIM(ymd)))
        READ(ymd,'(I2)') day0
        
        !
        ! Now extract from the time string
        !
        READ(hms(1:INDEX(hms,':')-1),'(I2)') hours0
        hms=TRIM(hms(INDEX(hms,':')+1:LEN_TRIM(hms)))
        READ(hms(1:INDEX(hms,':')-1),'(I2)') minutes0
        hms=TRIM(hms(INDEX(hms,':')+1:LEN_TRIM(hms)))
        READ(hms,'(I2)') seci
        !
        sec0 = hours0*3600. + minutes0*60. + seci
        CALL ymds2ju (year0, month0, day0, sec0, date0_tmp)
        date0_file(iff,tcnt) = date0_tmp
        !
        ! Now get the actual dates
        !
        tncstart(1) = 1
        tnccount(1) = nbtime_perfile(iff)
        IF ( check ) WRITE(*,*) "forcing_time : number of values read : ", tnccount(1)
        iret = NF90_GET_VAR(force_id(iff), time_id(iff,tcnt), time_read, tncstart, tnccount)
        IF (iret /= NF90_NOERR) THEN
           CALL ipslerr (3,'forcing_time',"An error occured while reading time from the file."," "," ")
        ENDIF
        !
        ! Convert the variable time from seconds since to julian days
        ! When we have a daily time step in the forcing a special procedure is used.
        ! Note that the time resolution could be 86400s or 1d in SRF/METEO/ISIMIP2B/
        ! But the time unit is not correctly defined in some time variables in this file
        !
        dttmp=(time_read(2)-time_read(1))*convtosec(iff)
        !
        DO it=1,nbtime_perfile(iff)
           !!time_infiles(tstart+it) = date0_file(iff,tcnt) + time_read(it)*convtosec(iff)/one_day
           IF ( dttmp < one_day ) THEN
              time_infiles(tstart+it) = date0_file(iff,tcnt) + time_read(it)*convtosec(iff)/one_day
           ELSE IF ( dttmp .GE. one_day .AND. convtosec(iff) .EQ. 1. ) THEN
              ! In the case of daily forcing, if the time resolution is 86400 and time unit is seconds 
              ! the start time has to be 00UTC in Julian days.
              time_infiles(tstart+it) = date0_file(iff,tcnt) + INT(time_read(it))/one_day
           ELSE IF ( dttmp .GE. one_day .AND. convtosec(iff) .EQ. one_day ) THEN
              ! In the case of daily forcing, and if the time resolution is 1d and time unit is day
              time_infiles(tstart+it) = date0_file(iff,tcnt) + INT(time_read(it))
           ENDIF
        ENDDO
        if ( check ) WRITE(*,*) "File ", iff, "goes from ",  time_infiles(tstart+1), " to ", &
             time_infiles(tstart+nbtime_perfile(iff))
        !
        ! Estimate the bounds as this information is not yet in the forcing file.
        !
        date_int = (time_infiles(tstart+nbtime_perfile(iff)) - time_infiles(tstart+1))/(nbtime_perfile(iff)-1)
        forcing_tstep_ave = date_int*one_day
        !
        ! If this is the first file then simply keep the name of the time axis. Else search the same name 
        ! in what has already been read
        !
        ! tcnt: corresponds to nbtax. There are three in the current unstructred grids, which is not useful.
        ! in writing the forcing file, the first time axis appeared should be the one to be used by physical variables
        ! nbtmethods: 4 cellmethods
        ! time_axename(nbtax*nbtmethods)
        !
        ! here have a problem. There are iterations over forcing files and then nbtax here.
        ! if it is the first file and if it is the first time axis, 
        IF ( iff == 1 ) THEN
           itbase = (tcnt-1)*nbtmethods
           time_axename(itbase+1:itbase+4) = timevarname
           WRITE(numout,*) "forcing_time for the first file : in ", iff, " found variable ", timevarname, tcnt, iff, itbase
           WRITE(numout,*) "forcing time axname first file ", time_axename(itbase+1:itbase+4)
           time_cellmethod(itbase+1) = "reference"
           time_cellmethod(itbase+2) = "start"
           time_cellmethod(itbase+3) = "cent"
           time_cellmethod(itbase+4) = "end"
        ELSE
           !
           ! If this is not the first file then we need to find the correct axis to add to.
           ! All information have already been saved with the first file.
           !
           DO ittmp=1,nbtax
              itbasetmp=(ittmp-1)*nbtmethods
              IF ( time_axename(itbasetmp+1) == timevarname ) THEN
                 itbase = itbasetmp
              ENDIF
           ENDDO
        ENDIF
        !
        ! Verify that the calendars match between forcing file and grid file
        !
        iret = NF90_GET_ATT(force_id(iff), time_id(iff,tcnt), "calendar", forcing_calendar)
        IF (iret == NF90_NOERR) THEN
           IF ( used_calendar == 'gregorian' ) THEN
              IF ( INDEX(lowercase(forcing_calendar), "gregorian") < 1 .AND. &
                   & INDEX(lowercase(forcing_calendar), "standard") < 1 ) THEN
                 CALL ipslerr (3,'forcing_time',"The model is set to a gregorian calendar but ",&
                      &          " does not match the calendar of the forcing file. ",&
                      &          " Check variable CalendarType in config.card")
              ENDIF
           ELSE IF ( used_calendar == 'noleap' ) THEN
              IF ( INDEX(lowercase(forcing_calendar), "noleap") < 1  .AND. &
                   & INDEX(lowercase(forcing_calendar), "365") < 1 ) THEN
                 CALL ipslerr (3,'forcing_time',"The model is set to a noleap calendar but ",&
                      &          " does not match the calendar of the forcing file. ",&
                      &          " Check variable CalendarType in config.card")
              ENDIF
           ELSE IF ( used_calendar == '360d' ) THEN
              IF ( INDEX(lowercase(forcing_calendar), "360") < 1 ) THEN
                 CALL ipslerr (3,'forcing_time',"The model is set to a 360d calendar but ",&
                      &          " does not match the calendar of the forcing file. ",&
                      &          " Check variable CalendarType in config.card")
              ENDIF
           ENDIF
           WRITE(numout,*) ">> Forcing calendar is ",forcing_calendar, " and matches the one for the model ", used_calendar
           
        ENDIF
        !
        ! Keep for future usage the various information on the time axis we have just read. This time axis can
        ! be understood in 3 different ways and each variable might use a different cell method for this time
        ! axis.
        !
        ! time(:,(tcnt-1)*nbtmethods+1) : corresponds to the reference time axis as it has been read from the file
        ! time(:,(tcnt-1)*nbtmethods+2) : is the time axis with a cell method which place the value at the 
        !                                beginning of the time interval
        ! time(:,(tcnt-1)*nbtmethods+3) : is the time axis corresponding to variables placed at the center of the 
        !                                time interval
        ! time(:,(tcnt-1)*nbtmethods+4) : for variables put at the end of the time interval over which they aere 
        !                                for instance averaged.
        !
        ! In variable time_cellmethod we will write the type of cell method as descirbed above so that the selection
        ! of the right axis for each variable can be made automaticaly.
        !
        ! We also keep the name of the time axis read in preparation of file where we might have to read more than one
        ! time axis.
        !
        ! it: ranges from 1 to the totoal number of data in the 3 forcing files (if there are three)
        ! time: time(nb_forcing_steps, nbtax*nbtmethods), having total number of forcing data over 3 full forcing files
        ! if there are 3 forcing files
        DO it=tstart+1,tstart+nbtime_perfile(iff)
           !
           ! Reference time
           !
           time(it,itbase+1) = time_infiles(it)
           time_bounds(it,itbase+1,1) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+1,2) = time_infiles(it)+date_int/2.0
           !
           ! Start cell method
           time(it,itbase+2) = time_infiles(it)+date_int/2.0
           time_bounds(it,itbase+2,1) = time_infiles(it)
           time_bounds(it,itbase+2,2) = time_infiles(it)+date_int
           !
           ! Centered cell method
           time(it,itbase+3) = time_infiles(it)
           time_bounds(it,itbase+3,1) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+3,2) = time_infiles(it)+date_int/2.0
           !
           ! End cell method
           time(it,itbase+4) = time_infiles(it)-date_int/2.0
           time_bounds(it,itbase+4,1) = time_infiles(it)-date_int
           time_bounds(it,itbase+4,2) = time_infiles(it)
           !
        ENDDO
        !
        ! Keep the number of the file from which we read this time.
        !
        time_sourcefile(tstart+1:tstart+nbtime_perfile(iff))=iff
        !
        IF ( check ) WRITE(*,*) "forcing_time : finished file ", iff
        !
     ENDDO  ! end of loup over nbtax
     !
     ! Before moving to the next file advance the pointer in the time arrays.
     ! tstart: time index value. the first value corresponds to the first value of the first forcing file
     tstart=tstart+nbtime_perfile(iff)
     !
  ENDDO ! end of loup over forcing files
  !
  IF ( check ) WRITE(*,*) "forcing_time : All files have been processed"
  !
  ! Verify that the forcing comes in regular time intervals. If not, many of the 
  ! interpolation schemes will fail.
  ! This is only done on the first time axis ... is it enough ?
  !
  DO ittmp=1,nbtax
     itbase=(ittmp-1)*nbtmethods
     !
     date_int = (time(nb_forcing_steps,itbase+1) - time(1,itbase+1))/(nb_forcing_steps-1)
     forcing_tstep_ave = date_int*one_day
     !
     timeint(:) = 0
     DO it=1, nb_forcing_steps-1
        timeint(it) = time(it+1,itbase+1)-time(it,itbase+1)
     ENDDO
     !
     IF (  MAXVAL(timeint(1:nb_forcing_steps-1)) > date_int+0.1*date_int .OR.&
          & MINVAL(timeint(1:nb_forcing_steps-1)) < date_int-0.1*date_int) THEN
        WRITE(*,*) "The time steping of the forcing files does not seem to be regular on axis",time_axename(itbase+1),":"
        WRITE(*,*) "Average time step : ", date_int, "days = ", date_int*one_day, "sec."
        imax = MAXLOC(timeint(1:nb_forcing_steps-1))
        imin = MINLOC(timeint(1:nb_forcing_steps-1))
        WRITE(*,*) "Maximum time step : ", MAXVAL(timeint(1:nb_forcing_steps-1)), " at ", imax(1)
        WRITE(*,*) "Minimum time step : ", MINVAL(timeint(1:nb_forcing_steps-1)), " at ", imin(1)
        WRITE(*,*) "++++ Values around Maximum"
        DO it=MAX(imax(1)-5,1),MIN(imax(1)+5,nb_forcing_steps)
           WRITE(*,*) it, " from file ", time_sourcefile(it), " Value ", time(it,itbase+1)
           CALL forcing_printdate(time(it,itbase+1), "Time values.")
        ENDDO
        WRITE(*,*) "++++ Values around Minimum"
        DO it=MAX(imin(1)-5,1),MIN(imin(1)+5,nb_forcing_steps)
           WRITE(*,*) it, " from file ", time_sourcefile(it), " Value ", time(it,itbase+1)
           CALL forcing_printdate(time(it,itbase+1), "Time values.")
        ENDDO
        CALL ipslerr (3,'forcing_time', 'The time handling could be improved to allow the driver',&
             & "to cope with irregular forcing files."," ")
     ENDIF
  ENDDO
  !
  ! Print some test values
  !
  DO ittmp=1,nbtax
     itbase=(ittmp-1)*nbtmethods
     !
     WRITE(*,*) "Bounds for axis ",time_axename(itbase+1)," :"
     !
     CALL forcing_printdate(time_bounds(1,itbase+1,1), "Start time of first forcing interval.")
     CALL forcing_printdate(time_bounds(1,itbase+1,2), "End time of first forcing interval.")
     CALL forcing_printdate(time_bounds(nb_forcing_steps,itbase+1,1), "Start time of last forcing interval.")
     CALL forcing_printdate(time_bounds(nb_forcing_steps,itbase+1,2), "End time of last forcing interval.")
  ENDDO
  !
  ! Keep the very first date of the time axis for future use
  ! 
  forcingstartdate = time(1,1)
  !
  ! Clean-up
  !
  DEALLOCATE(timeint, time_read)
  !
END SUBROUTINE forcing_time

! --------------------------------------------------------------------
! INTEGER FUNCTION  ExistVariable():
!    This function tests if the variable exists in the file.
! --------------------------------------------------------------------
FUNCTION  ExistVariable(fileindex, varname) RESULT(exists) 
  IMPLICIT  NONE
  INTEGER(i_std), INTENT(in) :: fileindex
  CHARACTER(LEN=*), INTENT(in) :: varname
  LOGICAL :: exists
  !
  CHARACTER(LEN=80) :: name
  INTEGER(i_std) :: iret, iv
  !
  exists=.FALSE.
  !
  DO iv=1,nvars(fileindex)
     !
     iret = NF90_INQUIRE_VARIABLE(force_id(fileindex), iv, name=name)
     !
     IF ( INDEX(name, varname) > 0 ) THEN
        exists = .TRUE.
     ENDIF
  ENDDO
  
END FUNCTION ExistVariable
  
!!  =============================================================================================================================
!! SUBROUTINE: forcing_meanprecip
!!
!>\BRIEF      Moyenne temporelle des precipitations totales du forcage, par maille.
!!
!! DESCRIPTION: Balaie tout l'axe temps de tous les fichiers de forcage et rend, pour chaque
!!              point de terre du processus, la moyenne de (Rainf + Snowf) convertie en mm/jour.
!!
!!              POURQUOI : au premier pas de temps d'un demarrage a froid, precip_longterm vaut
!!              zero, donc dynamic_parameters rabat pipe_tune2 sur son plancher min_pipe_tune2.
!!              Or c'est ce jour-la que prescribe dimensionne les peuplements. Des que la memoire
!!              de pluie s'amorce, pipe_tune2 rejoint sa valeur climatique et LE MEME BOIS est
!!              relu avec un diametre divise par ~1,8 : Nmax triple, le RDI d'etablissement tombe
!!              de 0,700 a ~0,23 sans qu'un seul arbre ne meure (diagnostic du 2026-07-29).
!!              Initialiser precip_longterm avec cette moyenne supprime le transitoire.
!!
!!              Le resultat est une moyenne, PAS un cumul : precip_longterm est en mm/jour et
!!              c'est dynamic_parameters qui multiplie par 365. Pluie ET neige sont sommees,
!!              comme le fait stomate (precip = (precip_rain + precip_snow)*one_day/dt_sechiba).
!!
!! \n
!_ ==============================================================================================================================
  SUBROUTINE forcing_meanprecip(model_guess, meanprec)
    !
    IMPLICIT NONE
    !
    !! 0.1 Input variables
    CHARACTER(LEN=8), INTENT(in) :: model_guess               !! Type de grille devine par globgrd
    !! 0.2 Output variables
    REAL(r_std), INTENT(out)     :: meanprec(nbpoint_proc)    !! Precipitation moyenne (mm/jour)
    !
    !! 0.4 Local variables
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: sum_full, mean_full, mean_proc
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: rainf_c, snowf_c
    INTEGER(i_std) :: ifile, it0, nct, ntot
    CHARACTER(LEN=80) :: cellmethod
    LOGICAL :: got_any
    !_ ================================================================================================================================
    !
    ALLOCATE(mean_proc(nbpoint_proc,1))
    mean_proc(:,:) = zero
    got_any = .FALSE.
    !
    IF ( is_root_prc ) THEN
       ALLOCATE(sum_full(nbpoint_loc,1), mean_full(nbpoint_loc,1))
       ALLOCATE(rainf_c(nbpoint_loc,slab_size), snowf_c(nbpoint_loc,slab_size))
       sum_full(:,:) = zero
       ntot = 0
       !
       DO ifile = 1,nb_forcefile
          it0 = 1
          DO WHILE ( it0 <= nbtime_perfile(ifile) )
             nct = MIN(slab_size, nbtime_perfile(ifile) - it0 + 1)
             rainf_c(:,:) = zero
             snowf_c(:,:) = zero
             !
             ! Rainf, ou "precip" (deja total) pour les forcages a pas de temps journalier.
             IF ( ExistVariable(ifile, "Rainf") ) THEN
                CALL forcing_varforslab(ifile, "Rainf", it0, nct, 1, model_guess, .FALSE., rainf_c, cellmethod)
                got_any = .TRUE.
             ELSE IF ( ExistVariable(ifile, "precip") ) THEN
                CALL forcing_varforslab(ifile, "precip", it0, nct, 1, model_guess, .FALSE., rainf_c, cellmethod)
                got_any = .TRUE.
             ENDIF
             !
             ! Snowf n'est ajoute QUE s'il est present ET distinct : quand le forcage ne porte
             ! que "precip", le total est deja dans rainf_c et ajouter une reconstruction de
             ! neige compterait la meme eau deux fois.
             IF ( ExistVariable(ifile, "Snowf") .AND. ExistVariable(ifile, "Rainf") ) THEN
                CALL forcing_varforslab(ifile, "Snowf", it0, nct, 1, model_guess, .FALSE., snowf_c, cellmethod)
             ENDIF
             !
             sum_full(:,1) = sum_full(:,1) + SUM(rainf_c(:,1:nct) + snowf_c(:,1:nct), DIM=2)
             ntot = ntot + nct
             it0 = it0 + nct
          ENDDO
       ENDDO
       !
       ! kg/m2/s -> mm/jour. Si rien n'a pu etre lu on rend zero : slowproc_init retombera
       ! alors sur le scalaire PRECIP_LONGTERM_INIT plutot que d'initialiser a zero, ce qui
       ! reproduirait exactement le defaut que cette routine corrige.
       IF ( ntot > 0 .AND. got_any ) THEN
          mean_full(:,1) = sum_full(:,1)/REAL(ntot,r_std)*one_day
       ELSE
          mean_full(:,1) = zero
       ENDIF
       !
       WRITE(numout,*) "forcing_meanprecip : ", ntot, " pas de temps balayes sur ", nb_forcefile, " fichier(s)"
       IF ( ntot > 0 .AND. got_any ) THEN
          WRITE(numout,*) "forcing_meanprecip : pluie+neige moyenne (mm/an) min/moy/max = ", &
               MINVAL(mean_full(:,1))*365., SUM(mean_full(:,1))/REAL(nbpoint_loc,r_std)*365., &
               MAXVAL(mean_full(:,1))*365.
       ELSE
          WRITE(numout,*) "forcing_meanprecip : AUCUNE variable de precipitation lisible -> repli sur le scalaire"
       ENDIF
    ELSE
       ALLOCATE(mean_full(1,1))
       mean_full(:,:) = zero
    ENDIF
    !
    CALL scatter(mean_full, mean_proc)
    meanprec(:) = mean_proc(:,1)
    !
    IF ( ALLOCATED(sum_full) ) DEALLOCATE(sum_full)
    IF ( ALLOCATED(mean_full) ) DEALLOCATE(mean_full)
    IF ( ALLOCATED(mean_proc) ) DEALLOCATE(mean_proc)
    IF ( ALLOCATED(rainf_c) ) DEALLOCATE(rainf_c)
    IF ( ALLOCATED(snowf_c) ) DEALLOCATE(snowf_c)
    !
  END SUBROUTINE forcing_meanprecip

!!  =============================================================================================================================
!! SUBROUTINE: forcing_varforslab
!!
!>\BRIEF      
!!
!! DESCRIPTION:	This subroutine will read the named variable and put it in the right place in the 
!!              slab of data kept in the memory of the driver. 
!!
!! \n
!_ ==============================================================================================================================
SUBROUTINE forcing_varforslab(fileindex, varname, timestart, timecount, inslabpos, model_guess, is_spinupacc, data, cellmethod)
  !
  ! This subroutine will read the named variable and put it in the right place in the 
  ! slab of data kept in the memory of the driver.
  !
  INTEGER(i_std), INTENT(in) :: fileindex
  CHARACTER(LEN=*), INTENT(in) :: varname
  INTEGER(i_std), INTENT(in) :: timestart, timecount, inslabpos
  REAL(r_std), INTENT(inout) :: data(nbpoint_loc,slab_size)
  CHARACTER(LEN=*), INTENT(out) :: cellmethod
  CHARACTER(len=8), INTENT(in) :: model_guess
  LOGICAL, INTENT(in) :: is_spinupacc
  !
  ! Local variables
  !
  INTEGER(i_std) :: varid, windid, windndims, it, ig, iv
  INTEGER(i_std) :: iret, ndims
  INTEGER(i_std), DIMENSION(4) :: start, count
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: tmp_slab
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: tmp_slab2d
  CHARACTER(LEN=80) :: name
  CHARACTER(LEN=80) :: coordinates ! temporary variable for attribut coordinates
  LOGICAL :: windzero
  
  !
  ! Allocate the temporary data array if not already available
  !
  IF ( compressed ) THEN
     IF ( .NOT. ALLOCATED(tmp_slab) ) ALLOCATE(tmp_slab(ncdfcount,slab_size))
  ELSE
     IF ( .NOT. ALLOCATED(tmp_slab2d) ) ALLOCATE(tmp_slab2d(iim_glo,jjm_glo,slab_size))
  ENDIF
  !
  ! Reset the counters and flags to forget the past !
  !
  varid=-1
  windid=-1
  windzero=.FALSE.
  !
  ! Find the variable in the file
  !
  DO iv=1,nvars(fileindex)
     !
     iret = NF90_INQUIRE_VARIABLE(force_id(fileindex), iv, name=name, ndims=it)
     !
     ! varname: the requested variable defined in forcing_readslab_root
     ! name: the variable names existing in the file
     ! if varname = Tair and name = Tair_qc, then varid of tair_qc will be used, not that of tair
     ! this will introduce an error. Add one more condition to avoid this case.
     IF ( ( INDEX(name, varname) > 0 ) .AND. ( (LEN_TRIM(name) .EQ. LEN_TRIM(varname)) ) ) THEN
        varid = iv
        ndims = it
     ENDIF
     IF ( (INDEX(name, "Wind") > 0) .AND. (LEN_TRIM(name) == LEN_TRIM("Wind")) ) THEN
        windid = iv
        windndims = it
     ENDIF
     !
  ENDDO
  !
  ! Treat some special cases and catch errors
  !
  IF ( varid < 0 ) THEN
     !
     ! If we requested a wind component but did not find it, it could be that only the module is available.
     ! If that is the case, then we use the module (windid) for the U component and set V top zero.
     !
     IF ( INDEX(varname, "Wind_E") > 0 ) THEN
        varid = windid
        ndims = windndims
        windzero = .FALSE.
     ELSE IF ( INDEX(varname, "Wind_N") > 0 ) THEN
        windzero = .TRUE.
     ELSE
        CALL ipslerr (3,'forcing_varforslab',"Could not find variable",varname," in file.")
     ENDIF
  ENDIF
  !
  ! If there is some variable to be read then do it
  !
  IF ( .NOT. windzero ) THEN
     !
     ! Get the attributes needed for intepretating the data
     !
     ! First get the cell method used for this variable
     iret = NF90_GET_ATT(force_id(fileindex), varid, 'cell_methods', cellmethod)
     IF (iret /= NF90_NOERR) THEN
        ! If the attribute is not found then we set undefined
        cellmethod="undefined"
     ENDIF
     ! Read also attribute coordinates that might contain information about the name of the time axis.
     ! If the forcing file is produced by XIOS, this is the case.
     ! Reading of the time axis name can fail if the attribute coordinates contains something else than the time axis name.
     iret = NF90_GET_ATT(force_id(fileindex), varid, 'coordinates', coordinates)
     IF (iret == NF90_NOERR) THEN
        ! The attribute 'coordinates' was found. This is probably containg the name of the time axis. 
        ! We add this information in the begining of the celllmethod variable.
        cellmethod=TRIM(coordinates)//":"//TRIM(cellmethod)
     ENDIF
     !
     !
     ! Getsize of data to be read from the netCDF file
     !
     !
     IF ( compressed ) THEN
        ! in the classical compressed grids such as crujra halfdeg, there are 67420 land point, with x 360 and y 180
        ! the Tair(tstep, land), nav_lon(y, x), contfrac(y, x), land(land)
        !
        ! in the unstructured grids of dynamico, there are 16002 total pixels and 5500 land pixels. It is not compressed.
        ! Tair(time, 16002), contfrac(16002), nav_lon(16002)
        ! we have to modify the following code in the case of dynamico
        !
        IF ( ndims == 2 ) THEN
           start = (/ncdfstart,timestart,0,0/)
           count = (/ncdfcount,timecount,0,0/)
        ELSE IF ( ndims == 3 ) THEN
           start = (/ncdfstart,1,timestart,0/)
           count = (/ncdfcount,1,timecount,0/)
        ELSE
           CALL ipslerr (3,'forcing_varforslab',"Compressed variable : ",varname,&
                &        " does not have a compatible number of dimensions.")
        ENDIF
        !
        iret = NF90_GET_VAR(force_id(fileindex), varid, tmp_slab, start, count)
        IF (iret /= NF90_NOERR) THEN
           CALL ipslerr (3,'forcing_varforslab',"Could not read from file variable : ",varname," Compressed in the file.")
        ENDIF
        !
        ! Zoom into the data and put it in the right place in the slab of data.
        ! for the case of unstructured grids, we call forcing_reindex2dt(nbin, tin, slab_in, nbout, tout, slab_out, tstart, reindex)
        !orig
        IF ( model_guess .EQ. 'unstruct' .AND. .NOT. is_spinupacc) THEN
        !   CALL forcing_reindex(ncdfcount, timecount, tmp_slab, nbpoint_loc, slab_size, data, inslabpos, reindex_loc)
        ! origind(iim_glo, jjm_glo), with the values being the indices of land points in reference to 5500 land points
        ! lindex_loc: dimension of 5500 points, with the values being the indices of land points in reference to 16002 total points
           CALL forcing_reindex2dt_unstruct(ncdfcount, timecount, tmp_slab, nbpoint_loc, slab_size, data, inslabpos, lindex_loc, origind(1,:))
        ELSE IF ( model_guess .EQ. 'unstruct' .AND. is_spinupacc) THEN
           CALL forcing_reindex2dt_spinup(ncdfcount, timecount, tmp_slab, nbpoint_loc, slab_size, data, inslabpos)
        ELSE
           CALL forcing_reindex(ncdfcount, timecount, tmp_slab, nbpoint_loc, slab_size, data, inslabpos, reindex_loc)
        ENDIF
        
     ELSE
        !
        IF ( ndims == 3 ) THEN
           start = (/1,1,timestart,0/)
           count = (/iim_glo,jjm_glo,timecount,0/)
        ELSE IF (ndims == 4 ) THEN
           start = (/1,1,1,timestart/)
           count = (/iim_glo,jjm_glo,1,timecount/) 
       ELSE
           CALL ipslerr (3,'forcing_varforslab',"Full lat Lon variable : ",varname,&
                &        " does not have a compatible number of dimensions.")
        ENDIF
        !
        !firstly read in the data organized in lon x lat dimension
        iret = NF90_GET_VAR(force_id(fileindex), varid, tmp_slab2d, start, count)
        IF (iret /= NF90_NOERR) THEN
           WRITE(*,*) TRIM(NF90_STRERROR(iret))
           WRITE(*,*) "File =", fileindex, "Size =", SIZE(tmp_slab2d,DIM=1), SIZE(tmp_slab2d,DIM=2), SIZE(tmp_slab2d,DIM=3)
           WRITE(*,*) "Start :", start(1:3)
           WRITE(*,*) "Count :", count(1:3)
           CALL ipslerr (3,'forcing_varforslab',"Could not read from file variable : ",varname," Not compressed.")
        ENDIF
        !
        ! then zoom into the data and put it in the right place in the slab of data, that is, into 1D dimension in space
        !
        CALL forcing_reindex(iim_glo, jjm_glo, timecount, tmp_slab2d, nbpoint_loc, slab_size, data, inslabpos, reindex2d_loc)
     ENDIF
  ELSE
     cellmethod="time: instantaneous"
     DO it=0,timecount-1
        data(:,inslabpos+it) = zero
     ENDDO
  ENDIF
  !
END SUBROUTINE forcing_varforslab

!!  =============================================================================================================================
!! SUBROUTINE: forcing_attributetimeaxe
!!
!>\BRIEF  Find the time axis which corresponds to the variable at hand.    
!!
!! DESCRIPTION:	 We interpret the cell_method provided in the netCDF file so that
!!               we can determine how we need to understand the values we read.
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE forcing_attributetimeaxe(cellmethod, timeindex)
  !
  ! We will analyse the time axis of the cell method found in the NetCDF file in order to
  ! attribute the correct time axis to this variable. 
  !
  CHARACTER(LEN=*), INTENT(in) :: cellmethod
  INTEGER(i_std), INTENT(out)  :: timeindex
  !
  INTEGER(i_std) :: itax, timepos, pos, lentime, itbase, im
  CHARACTER(LEN=20) :: TARGET, tmpstr
  CHARACTER(LEN=80) :: method
  !
  ! Clean the string to delete spaces in front of ":" and "("
  !
  method = cellmethod
  DO WHILE ( INDEX(method," :") > 0 )
     pos = INDEX(method," :")
     method = method(1:pos-1)//method(pos+1:LEN_TRIM(method))
  ENDDO
  DO WHILE ( INDEX(method,"( ") > 0 )
     pos = INDEX(method,"( ")
     method = method(1:pos)//method(pos+2:LEN_TRIM(method))
  ENDDO
  !
  ! Go through all the time axes we have to find the right one based on the information
  ! read from attribute cell_methods or coordinates
  !
  timeindex=0
  DO itax=1,nbtax
     !
     itbase=(itax-1)*nbtmethods
     ! find the time axis name in the cell method
     TARGET = TRIM(time_axename(itbase+1))//":"
     timepos = INDEX(method,TRIM(TARGET))
     !
     ! If we found the time axis then we look for the method with a time position description
     ! which is expected to be between parenthesis. For instance : mean(end)
     !
     IF ( timepos > 0 ) THEN
        !
        lentime=LEN_TRIM(time_axename(itbase+1))
        tmpstr = method(lentime+2:LEN_TRIM(method))
        !
        ! If there is ":" then there is information for another axis which needs to be deleted
        !
        IF ( INDEX(tmpstr,":") > 0 ) THEN
           tmpstr = tmpstr(1:INDEX(tmpstr,":")-1)
        ENDIF
        !
        ! Now that we have found a time axis see if we have between parenthesis a position
        ! on that time avis.
        !
        ! Look for a "("  
        IF ( INDEX(tmpstr, "(") > 0 ) THEN
           DO im=1,nbtmethods
              TARGET = "("//TRIM(time_cellmethod(itbase+im))
              timepos = INDEX(tmpstr,TRIM(TARGET))
              IF ( timepos > 0 ) THEN
                 timeindex = itbase+im
              ENDIF
           ENDDO
           !
        ELSE 
           ! There is no "(" so we set the centered axis by defaut
           DO im=1,nbtmethods
              IF ( INDEX(time_cellmethod(itbase+im), "cent") > 0 ) THEN
                 timeindex = itbase+im
              ENDIF
           ENDDO
        ENDIF
        !
        ! The name of the time axis was found but no method could be identified
        !
        IF ( timeindex < 1 ) THEN
           CALL ipslerr (3,'forcing_attributetimeaxe',"Found a time axis name but could not identify method.", &
                "This should not happen !"," ")
        ENDIF
        !
     ELSE
        ! Continue in loop over nbtax
     ENDIF
  ENDDO


  !
  ! If still no corresponding time axis name has be found, then we use the first centered one.
  ! This is the case if none of the attributes cell_methods or coordinates were found or if they 
  ! didn't include the name of the corresponding time axis.
  IF ( timeindex < 1 ) THEN
     CALL ipslerr(2,'forcing_attributetimeaxe',&
          'No information about the name of the time axis was found in the attributes cell_method or coordinantes', &
          'The first centered time axis will be used as defaut.', &
          'Attribute coordinates is optional but must not contain something else than the name of the time axis.')
     itax=1
     DO WHILE ( timeindex < 1 ) 
        IF ( INDEX(time_cellmethod(itax), "cent") > 0 ) THEN
           timeindex = itax
        ELSE
           itax = itax + 1
        ENDIF
     ENDDO
  END IF
  
END SUBROUTINE forcing_attributetimeaxe

!!  =============================================================================================================================
!! SUBROUTINE: forcing_filenamecheck
!!
!>\BRIEF   Check the list of files obtained from the calling program.   
!!
!! DESCRIPTION:	A small code to check the forcing files. They have to be NetCDF (i.e. .nc termination) and 
!!              we dont keep files that appear more than once in the list.  
!!
!! \n
!_ ==============================================================================================================================

SUBROUTINE forcing_filenamecheck(filenames_in, nb_files)
  !
  ! A small code to check the forcing files. They have to
  ! be NetCDF (i.e. .nc termination) and we dont keep files 
  ! that appear more than once in the list.
  !
  !
  ! INPUT
  !
  CHARACTER(LEN=*), DIMENSION(:), INTENT(in) :: filenames_in
  INTEGER(i_std), INTENT(out)                :: nb_files
  !
  ! LOCAL
  !
  INTEGER(i_std) :: ii, is, sizein
  LOGICAL        :: notfound
  !
  sizein = SIZE(filenames_in)
  IF ( sizein > 0 ) THEN
     IF ( ALLOCATED(forfilename) ) THEN
        DEALLOCATE(forfilename)
     ENDIF
     ALLOCATE(forfilename(sizein))
     nb_files=0
  ELSE
     CALL ipslerr (3,'forcing_filenamecheck',"List of forcing files is empty.","Please check your run.def file."," ")
  ENDIF
  !
  DO ii=1,sizein
     IF ( INDEX(filenames_in(ii), '.nc') > 0 ) THEN
        IF ( nb_files == 0 ) THEN
           nb_files = nb_files+1
           forfilename(nb_files)=TRIM(ADJUSTL(filenames_in(ii)))
        ELSE
           notfound=.TRUE.
           DO is=1,nb_files
              IF ( INDEX(TRIM(filenames_in(ii)), TRIM(ADJUSTL(forfilename(is)))) > 0 ) notfound=.FALSE.
           ENDDO
           IF ( notfound ) THEN
              nb_files = nb_files+1
              forfilename(nb_files)=TRIM(adjustl(filenames_in(ii)))
           ENDIF
        ENDIF
     ELSE
        !!! This is not a NetCDF file, so we ignore it
     ENDIF
  ENDDO
  !
  !
END SUBROUTINE forcing_filenamecheck

!!  =============================================================================================================================
!! FUNCTION: lowercase, FindMinimum, Swap 
!!
!>\BRIEF      Help functions for the forcing_tools module.
!!
!! DESCRIPTION:	  
!!
!! \n
!_ ==============================================================================================================================
function strlowcase(instr) result(outstr)
    ! turns a string to lower case
    character(len=*), intent(in)  :: instr
    character(len=LEN(instr))    :: outstr

    integer :: i, ascii, la, ua

    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= ua .and. ascii < ua+26) then ! in [A,Z]
          outstr(i:i) = achar( ascii - ua + la )
       endif
    enddo

    return
  end function strlowcase


FUNCTION lowercase(strIn) RESULT(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     IMPLICIT NONE

     CHARACTER(len=*), INTENT(in) :: strIn
     CHARACTER(len=LEN(strIn)) :: strOut
     INTEGER :: i,j

     DO i = 1, LEN(strIn)
          j = IACHAR(strIn(i:i))
          IF (j>= IACHAR("A") .AND. j<=IACHAR("Z") ) THEN
               strOut(i:i) = ACHAR(IACHAR(strIn(i:i))+32)
          ELSE
               strOut(i:i) = strIn(i:i)
          END IF
     END DO

END FUNCTION lowercase  
!
! Some help function found on Internet : http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
!
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------
INTEGER FUNCTION  FindMinimum(x, Start, END)
  IMPLICIT  NONE
  REAL(r_std), DIMENSION(1:), INTENT(IN) :: x
  INTEGER(i_std), INTENT(IN)             :: Start, END
  REAL(r_std)                            :: Minimum
  INTEGER(i_std)                         :: Location
  INTEGER(i_std)                         :: i

  Minimum  = x(Start)		! assume the first is the min
  Location = Start			! record its position
  DO i = Start+1, END		! start with next elements
     IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
        Minimum  = x(i)		!      Yes, a new minimum found
        Location = i            !      record its position
     ENDIF
  END DO
  FindMinimum = Location        	! return the position
END FUNCTION  FindMinimum
! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------
SUBROUTINE  Swap(a, b)
  IMPLICIT  NONE
  REAL(r_std), INTENT(INOUT) :: a, b
  REAL(r_std)                :: Temp

  Temp = a
  a    = b
  b    = Temp
END SUBROUTINE  Swap
! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------
SUBROUTINE  Sort(x, Size)
  IMPLICIT  NONE
  REAL(r_std), DIMENSION(1:), INTENT(INOUT) :: x
  INTEGER(i_std), INTENT(IN)                :: Size
  INTEGER(i_std)                            :: i
  INTEGER(i_std)                            :: Location

  DO i = 1, Size-1			! except for the last
     Location = FindMinimum(x, i, Size)	! find min from this to last
     CALL  Swap(x(i), x(Location))	! swap this and the minimum
  END DO
END SUBROUTINE  Sort

END MODULE forcing_tools
