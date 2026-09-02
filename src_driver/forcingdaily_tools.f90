!  ==============================================================================================================================\n
!  MODULE forcingdaily_tools : The general idea of this module is to re-generate a diurnal cycle of forcing variables based
!                              on the daily mean values and Tairmin and Tairmax. The approach is to generate a temporal sub-domain
!                              of size "szsbd" which will correspond to "nbdays" days. Thus "nbdays" daily means will be used
!                              to regenerate the diurnal cycles. Doing more than one day allows to use higher order interpolations
!                              in case it is needed and avoid discontinuities. This process is performed by forcingdaily_gensubd when
!                              ever we come to a new day. The the subroutine forcingdaily_getvalues will extract from the "nbdays" of
!                              reconstructed diurnal cycle the values ORCHIDEE needs.
!                              For most variables we have a specific subroutine to re-generate the diurnal cycle with some specific
!                              parameters which allow to adjust the process.
!
!               forcingdaily_gensubd :
!               forcingdaily_getvalues :
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
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_driver/forcingdaily_tools.f90 $
!! $Date: 2026-02-20 11:06:31 +0100 (ven. 20 févr. 2026) $
!! $Revision: 9385 $
!_ ================================================================================================================================
!!
MODULE forcingdaily_tools
  !
  USE defprec
  USE netcdf
  !
  USE ioipsl
  USE constantes
  USE solar
  USE qsat_moisture
  !
  USE mod_orchidee_para
  
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: forcingdaily_gensubd, forcingdaily_getvalues
  PUBLIC :: choice_qair_interpol, qair_interpol
  !
  ! This PARAMETER essentially manages the memory usage of the module as it
  ! determines how much of the forcing will be uploaded from the netCDF file into
  ! memory.
  !
  INTEGER(i_std), SAVE                              :: current_day = -1
  REAL(r_std), PARAMETER                            :: dusk_angle = 0.01
  INTEGER(i_std), PARAMETER                         :: nbdays = 3
  INTEGER(i_std), PARAMETER                         :: spreadprec = 7200    !! Time over which the precipitation should be distributed (in sec.)
  REAL(r_std), PARAMETER                            :: convprec_temp = 20.0 !! Temperature above which all precipitation is supposed to be convective.
                                                                            !! i.e. rainfall occurs only over spreadprec. Below rainfall will last longer.
  INTEGER(i_std), PARAMETER                         :: tmaxshift = 10800    !! How long after the solar noon should Tairmax occus ? Time in seconds. 
  INTEGER(i_std), SAVE                              :: seed = 7865439    
  !
  INTEGER(i_std), SAVE                              :: szsubd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)      :: time_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: tair_subd, qair_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: hurs_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: ztq_subd, zuv_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: rainf_subd, snowf_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: solarang_subd, swdown_subd, lwdown_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)    :: u_subd, v_subd, ps_subd
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:)  :: sinangles
  INTEGER(i_std), SAVE, ALLOCATABLE, DIMENSION(:,:) :: idusk, irise, inoon
  !
  LOGICAL, SAVE :: choice_qair_interpol = .TRUE.
  LOGICAL, SAVE :: qair_interpol = .TRUE.
CONTAINS
!!
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_getval
!!
!!
!>\BRIEF        Extracts the forcing values needed by ORCHIDEE from teh re-generated diurnal cycles.
!!
!! DESCRIPTION: As we are in a re-geneted diurnal cycle case, not a lot of precaution is taken to identify the value to be
!!              used for the next integration interval of ORCHIDEE. Simply the values located closest to the middle of the
!!              integration interval. This could be improved based on what is coded in forcing_tools.f90.     
!!  
!! \n
!_ ==============================================================================================================================
!  
  SUBROUTINE forcingdaily_getvalues(time_int, dt, ztq, zuv, tair, qair, rainf, snowf, &
       &                            swdown, lwdown, solarang, u, v, ps)
    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)  :: time_int(2)                            !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)  :: dt                                     !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    REAL(r_std), INTENT(out) :: ztq(:), zuv(:)
    REAL(r_std), INTENT(out) :: tair(:), qair(:), rainf(:), snowf(:)
    REAL(r_std), INTENT(out) :: swdown(:), lwdown(:), solarang(:)
    REAL(r_std), INTENT(out) :: u(:), v(:), ps(:)
    !
    ! LOCAL
    !
    INTEGER(i_std) :: imin(1), i, nbpt
    REAL(r_std)    :: tloc
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: qair_sat, hurs
    !
    tloc = (time_int(1)+time_int(2))/2.0
    imin = MINLOC(ABS((tloc - (time_subd(:)+current_day))))
    !
    nbpt = SIZE(qair)
    ALLOCATE(qair_sat(nbpt))
    ALLOCATE(hurs(nbpt))
    !
    !
    tair(:) = tair_subd(:,imin(1))
    hurs(:) = hurs_subd(:,imin(1))
    rainf(:) = rainf_subd(:,imin(1))
    snowf(:) = snowf_subd(:,imin(1))
    swdown(:) = swdown_subd(:,imin(1))
    lwdown(:) = lwdown_subd(:,imin(1))
    solarang(:) = solarang_subd(:,imin(1))
    ztq(:) = ztq_subd(:,imin(1))
    zuv(:) = zuv_subd(:,imin(1))
    u(:) = u_subd(:,imin(1))
    v(:) = v_subd(:,imin(1))
    ps(:) = ps_subd(:,imin(1))
    !
    ! Different options exist in order to get air humidity :
    ! - A linear interpolation of the specific humidity (not ideal !!).
    ! - A linear interpolation of the relative humidity which is then applied to the saturated
    !   value at air temperature (This should be preferred).
    !
    IF ( qair_interpol ) THEN
       qair(:) = qair_subd(:,imin(1))
    ELSE
       CALL qsatcalc (nbpt, tair, ps/100.0, qair_sat)
       IF ( MAXVAL(hurs) > un ) THEN
          ! hurs is certainly in %
          qair(:) = qair_sat(:)*hurs(:)/100.0
       ELSE
          ! Here hurs is a ratio.
          qair(:) = qair_sat(:)*hurs(:)
       ENDIF
    ENDIF
  END SUBROUTINE forcingdaily_getvalues
!!
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_gensubd
!!
!!
!>\BRIEF       generates the sub-diurnal cycle for a number of days around the current time step of the ORCHIDEE simulation. 
!!
!! DESCRIPTION: This routine only works when we start a new day so that not too much work is done. At each new day first the
!!              the diurnal evolution of the solar angle is computed and then from there all the rest is derived. For swdown
!!              it is a trivial process but for the other variables more complex procedures are used. For those variables
!!              nothing could be invented the daily mean value is places et the center of the day and a linear interpolation is
!!              used.  
!!  
!! \n
!_ ==============================================================================================================================
!
  SUBROUTINE forcingdaily_gensubd(time_int, dt, iim, jjm, lon, lat, gindex_proc, &
       &                          szdom, szslab, time_slab, ztq, zuv, tair, tairmin, tairmax, &
       &                          qair, hurs, rainf, snowf, swdown, lwdown, u, v, ps)

    !
    ! ARGUMENTS
    !
    REAL(r_std), INTENT(in)    :: time_int(2)                          !! The time interval over which the forcing is needed.
    REAL(r_std), INTENT(in)    :: dt                                   !! timestep, i.e. distance in seconds between time_int(1) and time_int(2)
    INTEGER(i_std), INTENT(in) :: szdom, szslab
    INTEGER(i_std), INTENT(in) :: iim, jjm                             ! Size of 2D domain
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm)           ! Longitude and latitude
    INTEGER(i_std), INTENT(in) :: gindex_proc(szdom)
    REAL(r_std), INTENT(in)    :: time_slab(szslab)
    REAL(r_std), INTENT(in)    :: ztq(szdom,szslab), zuv(szdom,szslab)
    REAL(r_std), INTENT(in)    :: tair(szdom,szslab), tairmin(szdom,szslab), tairmax(szdom,szslab)
    REAL(r_std), INTENT(in)    :: qair(szdom,szslab), hurs(szdom,szslab)
    REAL(r_std), INTENT(in)    :: rainf(szdom,szslab), snowf(szdom,szslab)
    REAL(r_std), INTENT(in)    :: swdown(szdom,szslab), lwdown(szdom,szslab)
    REAL(r_std), INTENT(in)    :: u(szdom,szslab), v(szdom,szslab), ps(szdom,szslab)
    !
    ! LOCAL
    !
    REAL(r_std)       :: tloc
    INTEGER(i_std)    :: it, i, ist, imin(1), imax(1), tmin(1), iday, j
    INTEGER(i_std)    :: stpday, half_subd
    REAL(r_std)       :: julian
    !
    ! Get the options chosen by the user
    !
    !Config Key   = DAILY_QAIR_INTERPOL
    !Config Desc  = Decide if qair from the daily forcing should be interpolated 
    !Config Def   = false
    !Config If    = 
    !Config Help  = The daily forcing can generate sub-diurnal qair either from the
    !Config         daily mean qair through an interpolation, or generate it from the
    !Config         relative humidity and the reconstructed tair. The second option is
    !Config         preferred so the default value for DAILY_QIAR_INTERPOL=false. 
    !Config Units = [-]
    !
    !
    IF ( choice_qair_interpol ) THEN
       qair_interpol = .TRUE.
       CALL getin_p('DAILY_QAIR_INTERPOL', qair_interpol)
    ELSE
       WRITE(*,*) "qair_interpol = ", qair_interpol
       CALL ipslerr(2, 'forcingdaily_gensubd', 'The forcing file has constrained our choices and the user', &
            &          'choice is not being used.', '')
    ENDIF
    !
    ! Set date to middle of requested interval.
    !
    tloc = (time_int(1)+time_int(2))/2.0
    !
    IF (INT(tloc) .NE. current_day) THEN
       !
       ! Save the date on which the 3 days are centered
       !
       current_day = INT(tloc)
       !
       ! The sub-diurnal cycle needs to be generated for the 3 days around the current time step.
       ! Allocate memory needed.
       !
       IF ( .NOT. ALLOCATED(time_subd) ) THEN
          !
          szsubd=INT(nbdays*one_day/dt)
          CALL random_seed()
          !
          ALLOCATE(time_subd(szsubd))
          ALLOCATE(tair_subd(szdom,szsubd))
          ALLOCATE(qair_subd(szdom,szsubd))
          ALLOCATE(hurs_subd(szdom,szsubd))
          ALLOCATE(rainf_subd(szdom,szsubd))
          ALLOCATE(snowf_subd(szdom,szsubd))
          ALLOCATE(swdown_subd(szdom,szsubd))
          ALLOCATE(lwdown_subd(szdom,szsubd))
          ALLOCATE(solarang_subd(szdom,szsubd))
          ALLOCATE(ztq_subd(szdom,szsubd))
          ALLOCATE(zuv_subd(szdom,szsubd))
          ALLOCATE(u_subd(szdom,szsubd))
          ALLOCATE(v_subd(szdom,szsubd))
          ALLOCATE(ps_subd(szdom,szsubd))
          ALLOCATE(sinangles(iim,jjm,szsubd))
          ALLOCATE(idusk(szdom,nbdays))
          ALLOCATE(irise(szdom,nbdays))
          ALLOCATE(inoon(szdom,nbdays))
       ENDIF
       !
       ! Number of sub-diurnal time steps per day.
       ! half_subd : number of days on either side of current date to be added.
       stpday = NINT(REAL(szsubd/nbdays))
       half_subd=INT(nbdays/2.0)
       !
       ! Generate time axis for the sub-diurnal domain, over the number of days we have (nbdays)
       !
       DO it=1,szsubd
          ! Time will be relative to the day at which we are being called.
          time_subd(it) = ((-half_subd*one_day)+(it-0.5)*dt)/one_day
       ENDDO
       !
       ! The following routines forcingdaily_xxx interpolate the forcing to
       ! a subdaily (30-minute) time axis by making use of the solar angles
       ! and radiation or in the case of precipitation a random distribution
       !
       ! swdown_subd
       CALL forcingdaily_solar(tloc, iim, jjm, lon, lat, gindex_proc, &
            &                       szdom, szslab, time_slab, swdown, &
            &                       szsubd, swdown_subd, sinangles, irise, inoon, idusk)
       !
       ! Temperature
       ! 
       CALL forcingdaily_tair(tloc, dt, szdom, nbdays, szslab, time_slab, tairmin, tairmax, &
            &                 irise, inoon, szsubd, tair_subd)
       !
       ! LWdown
       !
       CALL forcingdaily_lwdown(tloc, dt, szdom, szslab, time_slab, &
            &                   lwdown, szsubd, stpday, lwdown_subd)
       !
       ! Precipitation
       !
       CALL forcingdaily_precip(tloc, dt, iim, jjm, lon, lat, gindex_proc, &
            &                       szdom, szslab, time_slab, rainf, snowf, tair, &
            &                       szsubd, rainf_subd, snowf_subd)
       !
       ! For Qair, U, V, PS there is no clever way to do the interpolation so
       ! a simple linear interpolation is used. ztq and ztu are the heights of
       ! the variables t&q and u&v. hurs can be used when qair is not available
       ! from the forcing files.
       !
       qair_subd(:,:) = undef_sechiba
       hurs_subd(:,:) = undef_sechiba
       ztq_subd(:,:) = undef_sechiba
       zuv_subd(:,:) = undef_sechiba
       u_subd(:,:) = undef_sechiba
       v_subd(:,:) = undef_sechiba
       ps_subd(:,:) = undef_sechiba
       !
       ! xxx_subd is a variable for nbdays, where 48 locations
       ! are available for each day. The nbday values from the
       ! forcing file are placed in the middle of each day. Thus
       ! location 24, 72, and 120.
       !
       DO i=1,szdom
          DO iday=1,nbdays
             ist = (iday-1)*stpday+stpday/2
             julian = INT(time_subd(ist) + current_day)+0.5
             imin = MINLOC(ABS((julian - time_slab(1:szslab))))
             qair_subd(i,ist) = qair(i,imin(1))
             hurs_subd(i,ist) = hurs(i,imin(1))
             ztq_subd(i,ist) = ztq(i,imin(1))
             zuv_subd(i,ist) = zuv(i,imin(1))
             u_subd(i,ist) = u(i,imin(1))
             v_subd(i,ist) = v(i,imin(1))
             ps_subd(i,ist) = ps(i,imin(1))
          ENDDO
       ENDDO
       !
       ! Now the forcing has been corretly placed into a matrix with
       ! a time dimension of nbdays*48. The values that need to be
       ! interpolated have the value undef_sechiba. Linear interpolation
       ! will be used.
       !
       IF (qair_interpol) THEN
          ! If qair is to be interpolated, hurs is not available in the forcing.
          ! no need to interpolate hurs
          CALL forcingdaily_linint(szdom, szsubd, nbdays, qair_subd)
       ELSE
          ! If qair is not available hurs need to interpolated
          CALL forcingdaily_linint(szdom, szsubd, nbdays, hurs_subd)
       ENDIF
       !
       ! interpolate the remaining variables
       CALL forcingdaily_linint(szdom, szsubd, nbdays, ztq_subd)
       CALL forcingdaily_linint(szdom, szsubd, nbdays, zuv_subd)
       CALL forcingdaily_linint(szdom, szsubd, nbdays, u_subd)
       CALL forcingdaily_linint(szdom, szsubd, nbdays, v_subd)
       CALL forcingdaily_linint(szdom, szsubd, nbdays, ps_subd)
    ELSE
       ! Nothing to do as the sub-diurnal cycle has already been generated
    ENDIF
  END SUBROUTINE forcingdaily_gensubd
  
!!
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_linint
!!
!>\BRIEF        Does a linear interpolation of the variables already placed in the variable.
!!
!! DESCRIPTION: It will replace the undefined values (x >= undef_sechiba) by the interpolated value. Before entering this
!!              subroutine all time series need to be set to undef_sechiba with only the point between the interpolation should
!!              be applied filled with actual values.
!!  
!! \n
!_ ==============================================================================================================================
!
  SUBROUTINE forcingdaily_linint(xsz, tsz, nbdays, x)
    !!
    !! Arguments
    !!
    INTEGER(i_std), INTENT(in)    :: xsz, tsz, nbdays
    REAL(r_std), INTENT(inout)    :: x(xsz,tsz)
    !!
    !! Local
    !!
    INTEGER(i_std) :: ii, it, i
    INTEGER(i_std) :: find, nind, sind
    REAL(r_std)    :: fval, nval, del
    !!
    ! xsz: the number of pixels
    DO ii=1,xsz
       !
       find = -1
       fval = undef_sechiba
       nind = 1
       nval = undef_sechiba
       !
       ! tsz: 144 time steps, 3 days
       ! nind: time step index
       ! x:input physical var to be interpolated
       DO WHILE (nind < tsz)
          !
          ! Look for the next defined value
          !
          DO WHILE (nval >= undef_sechiba .AND. nind < tsz)
             IF (x(ii,nind) < undef_sechiba) THEN
                nval=x(ii,nind)
             ELSE
                nind=nind+1
             ENDIF
             IF (nind >= tsz) THEN
                nval = undef_sechiba
                nind = tsz
             ENDIF
          ENDDO
          !
          ! Do the filling or interpolation between find and nind
          !
          IF ( find < 0 ) THEN
             DO i=1,nind
                x(ii,i) = nval
             ENDDO
          ELSE IF (nind == tsz) THEN
             DO i=find,tsz
                x(ii,i) = fval
             ENDDO
          ELSE
             del = (nval-fval)/(nind-find)
             DO i=find,nind
                x(ii,i) = fval+del*(i-find)
             ENDDO
          ENDIF
          
          !
          ! Move information to first index
          !
          find = nind
          fval = nval
          nval = undef_sechiba
          nind = nind+1
       ENDDO
    ENDDO
  END SUBROUTINE forcingdaily_linint
!!
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_solar
!!
!!
!>\BRIEF        Computes the diurnal cycle of the solarangle and incident solar radiation.
!!
!! DESCRIPTION: This is very close to what is done in forcing_tools.f90 for the SWdown interpolation. The added code here
!!              is to compute the time indicis of sun rise, noon and sun set. These will be important point in the
!!              re-generation of the diurnal cycles.
!!  
!! \n
!_ ==============================================================================================================================
!  
  SUBROUTINE forcingdaily_solar(tloc, iim, jjm, lon, lat, gindex_proc, &
       &                        szdom, szslab, time_slab, swdown_loc, &
       &                        szsubd, swdown_subd, sinangles, irise, inoon, idusk)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)    :: tloc
    INTEGER(i_std), INTENT(in) :: szdom, szslab
    INTEGER(i_std), INTENT(in) :: iim, jjm                             ! Size of 2D domain
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm)           ! Longitude and latitude
    INTEGER(i_std), INTENT(in) :: gindex_proc(szdom)
    REAL(r_std), INTENT(in)    :: time_slab(szslab)
    REAL(r_std), INTENT(in)    :: swdown_loc(szdom,szslab)
    INTEGER(i_std), INTENT(in) :: szsubd
    REAL(r_std), INTENT(out)   :: swdown_subd(szdom,szsubd)
    REAL(r_std), INTENT(out)   :: sinangles(iim,jjm,szsubd)
    INTEGER(i_std), INTENT(out)   :: idusk(szdom,nbdays)
    INTEGER(i_std), INTENT(out)   :: irise(szdom,nbdays)
    INTEGER(i_std), INTENT(out)   :: inoon(szdom,nbdays)    
    !
    ! Local
    !
    REAL(r_std)       :: cval, lval
    INTEGER(i_std)    :: it, ii, jj, i, ist, imin(1), imax(1), tmin(1), iday
    INTEGER(i_std)    :: year, month, day, hours, minutes
    INTEGER(i_std)    :: stpday, half_subd
    REAL(r_std)       :: sec, julian
    REAL(r_std), SAVE :: solaryearstart
    REAL(r_std)       :: sinang(iim,jjm), mean_sinang(iim,jjm,nbdays)
    INTEGER(i_std)    :: nbval(nbdays)
    LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: mask
    !
    ! Allocate memory
    !
    IF ( .NOT. ALLOCATED(mask) ) THEN
       ALLOCATE(mask(szsubd))
    ENDIF
    !
    ! Get some basic dates and dimensions
    !
    stpday = NINT(REAL(szsubd/nbdays))
    half_subd=INT(nbdays/2.0)
    CALL ju2ymds (tloc, year, month, day, sec)
    CALL ymds2ju (year, 1, 1, 0.0, solaryearstart)
    !
    mean_sinang(:,:,:) = 0.0
    nbval(:) = 0
    !
    ! Compute all solar angles
    !
    DO it=1,szsubd
       julian=time_subd(it) + current_day
       iday = INT(time_subd(it)+(nbdays-half_subd))
       CALL solarang (julian, solaryearstart, iim, jjm, lon, lat, sinang)
       DO ii=1,iim
          DO jj=1,jjm
             IF ( sinang(ii,jj) > zero .AND.  sinang(ii,jj) < dusk_angle ) THEN
                sinang(ii,jj) = dusk_angle
             ENDIF
             mean_sinang(ii,jj,iday) = mean_sinang(ii,jj,iday)+sinang(ii,jj)
          ENDDO
       ENDDO
       nbval(iday) = nbval(iday)+1
       !
       sinangles(:,:,it) = sinang(:,:)
       !
    ENDDO
    DO it=1,nbdays
       mean_sinang(:,:,it) = mean_sinang(:,:,it)/nbval(it)
    ENDDO
    !
    tmin = MINLOC(ABS((tloc-current_day)-time_subd(1:szsubd)))
    imin = MINLOC(ABS((tloc - time_slab(1:szslab))))
    sinang(:,:) = sinangles(:,:,tmin(1))
    !
    ! Set default values for irise, inoon and idusk so that in polar
    ! regions we are not left with these indicis outside of range.
    ! These values only matter when the sun never rises (sinangles < dusk_angle).
    !
    DO iday=1,nbdays
       irise(:,iday) = (iday-1)*stpday + INT(stpday/3)
       inoon(:,iday) = (iday-1)*stpday + INT(stpday/2)
       idusk(:,iday) = (iday-1)*stpday + INT(stpday/2) + INT(stpday/3)
    ENDDO
    !
    DO it=1,szsubd
       DO i=1,szdom
          !
          iday = INT(time_subd(it)+(nbdays-half_subd))
          ! Put Julian date to mid-day of current day
          julian = INT(time_subd(it) + current_day)+0.5
          imin = MINLOC(ABS((julian - time_slab(1:szslab))))
          !
          jj = ((gindex_proc(i)-1)/iim)+1
          ii = (gindex_proc(i)-(jj-1)*iim)
          !
          IF ( mean_sinang(ii,jj,iday) > zero ) THEN
             swdown_subd(i,it) = swdown_loc(i,imin(1))*sinangles(ii,jj,it)/mean_sinang(ii,jj,iday)
          ELSE
             swdown_subd(i,it) = zero
          ENDIF
          !
          !
          solarang_subd(i,it) = sinangles(ii,jj,it)
          !
          lval = sinangles(ii,jj,MAX(it-1,1))
          cval = sinangles(ii,jj,it)
          IF ( lval .LE.  dusk_angle .AND. cval .GT. dusk_angle ) THEN
             irise(i,iday) = it
          ENDIF
          IF ( lval .GT. dusk_angle .AND. cval .LE. dusk_angle ) THEN
             idusk(i,iday) = it
          ENDIF
       ENDDO
    ENDDO
    !
    ! Position the solar noon in each day
    !
    DO i=1,szdom

       jj = ((gindex_proc(i)-1)/iim)+1
       ii = (gindex_proc(i)-(jj-1)*iim)
       
       DO it=1,nbdays
          mask=.FALSE.
          DO ist=1,stpday
             mask((it-1)*stpday+ist) = .TRUE.
          ENDDO
          imax = MAXLOC(swdown_subd(i,:), mask)
          IF ( sinangles(ii,jj,imax(1)) > dusk_angle ) THEN
             inoon(i,it) = imax(1)
          ENDIF
       ENDDO
    ENDDO
    !
  END SUBROUTINE forcingdaily_solar
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_lwdown
!!
!!
!>\BRIEF        Re-generates the LWdown variable.
!!
!! DESCRIPTION: For the moment LW down is maintained to the mean value over the full diurnal cycle. An idea would be to
!!              to introduce a diurnal cycle related to the lower atmospheric temperature (Tair).
!!  
!! \n
!_ ==============================================================================================================================
!
!! Interpolation of LWdown
!!
  SUBROUTINE forcingdaily_lwdown(tloc, dt, szdom, szslab, time_slab, &
       &                         lwdown, szsubd, stpday, lwdown_subd)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)    :: tloc, dt
    INTEGER(i_std), INTENT(in) :: szdom, szslab, stpday, szsubd
    REAL(r_std), INTENT(in)    :: time_slab(szslab)
    REAL(r_std), INTENT(in)    :: lwdown(szdom,szslab)
    REAL(r_std), INTENT(out)   :: lwdown_subd(szdom,szsubd)
    !
    ! Local
    !
    INTEGER(i_std) :: i, iday, it, ist, imin(1)
    REAL(r_std)    :: julian
    !
    lwdown_subd(:,:) = undef_sechiba
    DO i=1,szdom
       DO iday=1,nbdays
          DO it=1,stpday
             ist = (iday-1)*stpday+it
             julian = INT(time_subd(ist) + current_day)+0.5
             imin = MINLOC(ABS((julian - time_slab(1:szslab))))
             lwdown_subd(i,ist) = lwdown(i,imin(1))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE forcingdaily_lwdown
  !! =============================================================================================================================
!! SUBROUTINE: forcingdaily_tair
!!
!!
!>\BRIEF        Re-generates the Tair variable.
!!
!! DESCRIPTION: 
!!  
!! \n
!_ ==============================================================================================================================
!
!! Interpolation of Tair
!!
  SUBROUTINE forcingdaily_tair(tloc, dt, szdom, nbdays, szslab, time_slab, tairmin, tairmax, &
       &                       irise, inoon, szsubd, tair_subd)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)    :: tloc, dt
    INTEGER(i_std), INTENT(in) :: szdom, szslab, nbdays, szsubd
    REAL(r_std), INTENT(in)    :: time_slab(szslab)
    REAL(r_std), INTENT(in)    :: tairmin(szdom,szslab), tairmax(szdom,szslab)
    INTEGER(i_std), INTENT(in) :: irise(szdom,nbdays), inoon(szdom,nbdays)
    REAL(r_std), INTENT(out)    :: tair_subd(szdom,szsubd)
    !
    ! Local
    !
    INTEGER(i_std) :: i, it, imin(1), in, ir
    INTEGER(i_std) :: stpday, nbshift
    REAL(r_std)    :: julian
    !
    stpday = NINT(REAL(szsubd/nbdays))
    nbshift = NINT(REAL(tmaxshift/dt))
    !
    tair_subd(:,:) = undef_sechiba
    DO i=1,szdom
       DO it=1,nbdays
          !
          ! Tairmax position
          !
          in = (it-1)*stpday + MOD(inoon(i,it)+nbshift-1, stpday)+1
          !! in = inoon(i,it)+nbshift
          julian = INT(time_subd(in) + current_day)+0.5
          imin = MINLOC(ABS((julian - time_slab(1:szslab))))
          tair_subd(i,in) = tairmax(i,imin(1))
          !
          ! Tairmin position
          !
          ir = irise(i,it)
          julian = INT(time_subd(ir) + current_day)+0.5
          imin = MINLOC(ABS((julian - time_slab(1:szslab))))
          tair_subd(i,ir) = tairmin(i,imin(1))
       ENDDO
    ENDDO
    !
    CALL forcingdaily_linint(szdom, szsubd, nbdays, tair_subd)
    !
  END SUBROUTINE forcingdaily_tair
!!
!! =============================================================================================================================
!! SUBROUTINE: forcingdaily_precip
!!
!>\BRIEF        Distributes the rainfall of a day on a period spreadprec in a random place within the diurnal cycle.
!!
!! DESCRIPTION: Rainfall is distributed randomly within the day over the spreadprec period. A wide room for improvements here.
!!              We need to think about the geographical variations of spreadprec and the most lieky time of day for rainfall.
!!  
!! \n
!_ ==============================================================================================================================
!  
  SUBROUTINE forcingdaily_precip(tloc, dt, iim, jjm, lon, lat, gindex_proc, &
       &                        szdom, szslab, time_slab, rainf, snowf, tair, &
       &                        szsubd, rainf_subd, snowf_subd)
    !
    ! Arguments
    !
    REAL(r_std), INTENT(in)    :: tloc, dt
    INTEGER(i_std), INTENT(in) :: szdom, szslab
    INTEGER(i_std), INTENT(in) :: iim, jjm                             ! Size of 2D domain
    REAL(r_std), INTENT(in)    :: lon(iim,jjm), lat(iim,jjm)           ! Longitude and latitude
    INTEGER(i_std), INTENT(in) :: gindex_proc(szdom)
    REAL(r_std), INTENT(in)    :: time_slab(szslab)
    REAL(r_std), INTENT(in)    :: rainf(szdom,szslab), snowf(szdom,szslab)
    REAL(r_std), INTENT(in)    :: tair(szdom,szslab)
    INTEGER(i_std), INTENT(in) :: szsubd
    REAL(r_std), INTENT(out)   :: rainf_subd(szdom,szsubd), snowf_subd(szdom,szsubd)
    !
    ! Local
    !
    INTEGER(i_std)   :: i, ist, iday, it, ip
    INTEGER(i_std)   :: imin(1)
    INTEGER(i_std)   :: stpday, nbstep
    REAL(r_std)      :: julian, rr, rainlength
    ! modified 202506 
    INTEGER :: seed_size, iseed
    INTEGER, ALLOCATABLE :: seed(:)
    REAL(r_std), ALLOCATABLE :: tair_new(:,:)
    REAL(r_std) :: temperature
    
    CALL RANDOM_SEED(size=seed_size)
    IF (.NOT. ALLOCATED(seed)) ALLOCATE(seed(seed_size))
    IF (.NOT. ALLOCATED(tair_new)) ALLOCATE(tair_new(szdom, nbdays))
    !
    stpday = NINT(REAL(szsubd/nbdays))
    DO i=1,szdom
       DO iday=1,nbdays
       	  ist = (iday-1)*stpday+INT(stpday/2)
          julian = INT(time_subd(ist) + current_day)+0.5
          imin = MINLOC(ABS((julian - time_slab(1:szslab))))
          tair_new(i, iday)= tair(i, imin(1))
       ENDDO
    ENDDO

    ! define seed by using temperature, considering high entropy
    DO iseed = 1, seed_size
       temperature = tair_new(MOD(iseed, szdom) + 1, MOD(iseed, nbdays) + 1)
       seed(iseed) = MOD(INT(temperature * 1000.0) + iseed * 37, HUGE(0))
    END DO
    CALL RANDOM_SEED(put = seed)

    ! nbdays 3, szsubd 144 (3day*24*2)
    ! stpday: the number of time steps each day, 48
    !stpday = NINT(REAL(szsubd/nbdays))
    !
    rainf_subd(:,:) = zero
    snowf_subd(:,:) = zero
    !
    DO i=1,szdom
       DO iday=1,nbdays
          !
          ist = (iday-1)*stpday+INT(stpday/2)
          julian = INT(time_subd(ist) + current_day)+0.5
          imin = MINLOC(ABS((julian - time_slab(1:szslab))))
          !
          ! Air temperature decides how long the rainfall will be !
          ! When it is cold the precipitation will last 3 times longer.
          !
          rainlength=spreadprec + MIN(MAX(convprec_temp - (tair(i,imin(1)) - tp_00), 0.0), &
               convprec_temp)/convprec_temp*3.0*spreadprec
          nbstep = NINT(REAL(rainlength/dt))
          ! we can define the seed by using the rainlenth
          ! But for all temperatures above 20°C the same seed would be used
          ! seed = INT(rainlength)
          ! Or, we can define seed by using tair to avoid the above issue
          !seed(1) = INT(tair(i,imin(1))*1000)
          !CALL RANDOM_SEED(put=seed)
          !
          IF (rainf(i,imin(1)) > zero .OR. snowf(i,imin(1)) > zero) THEN
             CALL random_number(rr)
             it = INT(rr*stpday)+1
             DO ist=1,nbstep
                ip = MOD(it-1+(ist-1),stpday)+1
                rainf_subd(i,(iday-1)*stpday+ip) = rainf(i,imin(1))*one_day/rainlength*dt
                snowf_subd(i,(iday-1)*stpday+ip) = snowf(i,imin(1))*one_day/rainlength*dt
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    IF (ALLOCATED(tair_new)) DEALLOCATE(tair_new)
    IF (ALLOCATED(seed)) DEALLOCATE(seed)
    !
  END SUBROUTINE forcingdaily_precip
!
END MODULE forcingdaily_tools
