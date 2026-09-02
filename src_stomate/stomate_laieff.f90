! =================================================================================================================================
! MODULE       : stomate_laieff
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that are used to calculate the effective
!!             LAI.  Right now this is only used in the albedo routines.  I
!!             originally tried to put these routines there, but that creates
!!             some dependency problems that could not be resolved, since
!!             the effective LAI is calculated in stomate.f90 along with
!!             the LAI.
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_laieff.f90 $
!! $Date: 2025-03-05 19:40:39 +0100 (mer. 05 mars 2025) $
!! $Revision: 8845 $
!! \n
!_ ================================================================================================================================

MODULE stomate_laieff

  ! Modules used:
  USE xios_orchidee
  USE constantes
  USE pft_parameters
  USE structures
  USE ioipsl_para, ONLY : ipslerr_p
  USE stomate_stand_structure
  USE stomate_data
  USE function_library, ONLY: wood_to_height_eff, wood_to_height, &
                              biomass_to_lai, cc_to_lai, wood_to_cn_dia, &
                              get_printlev
  USE matrix_resolution

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC effective_lai, fitting_laieff, calculate_z_level_photo, &
       find_lai_per_level, combine_lai_levels, &
       calculate_laieff_fit, stomate_laieff_initialize, partial_spheroid_vol

  INTEGER(i_std), SAVE :: printlev_loc                !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
  LOGICAL, SAVE        :: firstcall_laieff = .TRUE.   !! First time a subroutine in current module is called
!$OMP THREADPRIVATE(firstcall_laieff)
  LOGICAL, SAVE        :: hack_pwc_h                  !! Hack to bypass an error on variabl pwc_h in effectiv_lai subroutine
!$OMP THREADPRIVATE(hack_pwc_h)


CONTAINS


!! ================================================================================================================================
!! SUBROUTINE     : stomate_laieff_initialize
!!
!>\BRIEF        : Initialize module stomate_laieff
!!
!! DESCRIPTION : This subroutine is called from stomate_initialize. It will initialize printlev_loc variable by 
!!               reading PRINTLEV_stomate_laieff from run.def.
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_laieff_initialize()
    
    !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    printlev_loc=printlev
    firstcall_laieff=.FALSE.
    
    !Config Key   = HACK_PWC_H
    !Config Desc  = Allow error on varible pwc_h and set it to 1
    !Config If    = 
    !Config Def   = n
    !Config Help  = 
    !Config Units = [y/n]   
    hack_pwc_h=.FALSE.
    CALL getin_p('HACK_PWC_H',hack_pwc_h)
    
  END SUBROUTINE stomate_laieff_initialize
    

!! ================================================================================================================================
!! SUBROUTINE 	: calculate_z_level_photo
!!
!>\BRIEF        
!!
!! DESCRIPTION  :  This is routine to calculate the physical heights of
!!                 all the levels used to calculate photosynthesis.  The bottom layer
!!                 will always have a height of zero, since it must be at the ground
!!                 level.  I tried to make all the layers evenly spaced starting from
!!                 zero to the top of the canopy, but this does not work very well
!!                 for tall trees; all the vegetation ends up in the top level and
!!                 nothing in the lower levels.  So now I try to find the bottom of
!!                 the canopy for trees as well, and start the second level from there.
!!                 The top level is always just below the top of the canopy (unless there
!!                 is only one level).  If there are multiple levels used in the energy
!!                 budget, the spacing between the levels will not be uniform.
!!
!!                 If there is no leaf biomass, all the levels are set to zero since
!!                 they should not be used and it would be a waste of time to calculate.
!!
!!                 Please be aware of the vertical layers. The original design
!!                 using nlevels and nlai here in
!!                 calculate_z_level_photo has only been tested together with
!!                 the new multi-layer energy budget, when set to nlevels=1 and
!!                 nlai=jnlvls_canopy+1.  
!!    NOTE:
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::z_level_photo
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE calculate_z_level_photo(npts, circ_class_biomass, circ_class_n, z_level_photo)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(IN)                        :: npts                 !! Domain size - number of pixels (unitless)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN)  &
                                                     :: circ_class_biomass   !! Biomass of the componets of the model  
                                                                             !! tree within a circumference
                                                                             !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)        :: circ_class_n         !! Number of trees within each circ class
                                                                             !! class @tex $(m^{-2})$ @endtex

    !! 0.2 Output variables
    REAL(r_std),DIMENSION(:,:,:),INTENT(OUT)         :: z_level_photo        !! the height of the bottom of each of our levels 
                                                                             !! @tex $(m)$ @endtex
                                                                             !! note that level 1 is closest to the ground


    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER :: ipts,ivm,il,il_photo,icir                                     !! Indices (unitless)
    REAL(r_std),DIMENSION(ncirc)  :: heights                                 !! Tree heights @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(ncirc)  :: cn_dia                                  !! Crown diameters @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(ncirc)  :: cn_bottom                               !! Height of the bottom of the crown @tex $(m)$ @endtex
    REAL(r_std) :: max_height                                                !! Highest actual vegetation height @tex $(m)$ @endtex
    REAL(r_std) :: min_height                                                !! Lowest bottom of the canopy @tex $(m)$ @endtex
    REAL(r_std) :: upper_level_boundary                                      !! The height of the top of the current energy budget
                                                                             !! level @tex $(m)$ @endtex
    REAL(r_std) :: lower_level_boundary                                      !! The height of the bottom of the current energy budget
                                                                             !! level @tex $(m)$ @endtex
    REAL(r_std) :: ideal_level_thickness                                     !! The thickness of the photosynthesis levels in this
                                                                             !! energy budget level if there are no other concerns.
                                                                             !! @tex $(m)$ @endtex
    REAL(r_std) :: block_thickness                                           !! The thickness of a chunk of PS levels together.
                                                                             !! @tex $(m)$ @endtex
    REAL(r_std) :: canopy_thickness                                          !! The thickness of the whole vegetation layer.
                                                                             !! @tex $(m)$ @endtex
    REAL(r_std) :: canopy_center                                             !! The height of the center of the whole vegetation layer.
                                                                             !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(ncirc) :: dummy                                    !! Non-optional output argument that is not used.   
!_ ================================================================================================================================

    IF (firstcall_laieff) CALL ipslerr_p (3,'calculate_z_level_photo', 'Initialization is not done before first call','','')

    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering calculate_z_level_photo'

    ! Debug
    IF (printlev_loc.GE.4) WRITE(numout,*) 'Entering calculate_z_level_photo',nlevels,nlevels_tot

    ! If we are not using any additional levels for 
    ! the photosynthesis, this is easy.
    IF(nlai .EQ. 1) THEN
       DO il=1,nlevels

          z_level_photo(:,:,il)=z_level(il)

          ! Debug
          IF(printlev_loc >= 4)THEN
             WRITE(numout,*)'is this were we calculate z_level_photo?'
             WRITE(numout,*) 'yes, z_level_photo is', z_level_photo(:,:,il)
          ENDIF
          !-

       ENDDO
       RETURN
    ENDIF

    IF(nlevels .GT. 1)THEN
       WRITE(numout,*) 'WARNING: calculate_z_level_photo has not been tested yet'
       WRITE(numout,*) 'WARNING: for the multilevel energy budget.  Be sure to '
       WRITE(numout,*) 'WARNING: do this and then remove this warning message.'
    ENDIF

    ! This is just a quick test.  The difference between the energy level heights
    ! should not be so small that it is impossible to fit all our photosynthesis
    ! levels in.
    DO il=1,nlevels-1

       ideal_level_thickness=(z_level(il+1)-z_level(il))&
            /REAL(nlai,r_std)

       ! Debug
       IF(printlev_loc.GT.4)THEN
          WRITE(numout,*)'ideal_level_thickness',ideal_level_thickness
       END IF
       !-

       DO ivm=1,nvm

          ! Debug
          IF(printlev_loc.GT.4)THEN
             WRITE(numout,*)'min_level_sep,ivm', min_level_sep(ivm),ivm
          END IF
          !-
   
          ! Error checking
          IF(ideal_level_thickness .LE. min_level_sep(ivm))THEN
             WRITE(numout,*) 'ERROR: We do not have enough space between the energy '
             WRITE(numout,*) '       levels to fit all the photosynthesis levels you '
             WRITE(numout,*) '       are asking for.  Either reduce nlai or min_level_sep, '
             WRITE(numout,*) '       or increase the width of the levels for the energy budget.'
             WRITE(numout,*) 'il, ivm: ',il,ivm
             WRITE(numout,*) 'ideal_level_thickness, min_level_sep(ivm): ',&
                  ideal_level_thickness, min_level_sep(ivm)
             WRITE(numout,*) 'z_level(il+1),z_level(il): ',z_level(il+1),z_level(il)
             WRITE(numout,*) 'nlai: ',nlai
             CALL ipslerr_p (3,'calculate_z_level_photo', 'Not enough space.','','')
          ENDIF
          !-

       ENDDO

    ENDDO

    DO ipts=1,npts
       DO ivm=1,nvm
          ! Skip all of this if we have no leaf biomass
          IF(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)) .LE. min_stomate)THEN
             z_level_photo(ipts,ivm,:)=zero
             CYCLE
          ENDIF

          IF(is_tree(ivm))THEN

             ! We need to know the maximum height of the vegetation so that
             ! we don't calculate a lot of levels with will remain empty.
             heights(:)=wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                  ivm,pipe_tune2(ipts,ivm))
             max_height=MAXVAL(heights(:))

             ! Find the vertical crown diameter and subtract this from the height,
             ! so that we know where the bottom of the canopy is.
             CALL wood_to_cn_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                  circ_class_n(ipts,ivm,:),ivm,dummy,cn_dia(:),&
                  pipe_tune2(ipts,ivm))
             cn_bottom(:)=heights(:)-cn_dia(:)
             min_height = MINVAL(cn_bottom(:))
             IF(min_height < min_stomate) min_height = 0.0001_r_std

             ! Notice that after coppicing, we may have zero aboveground wood, so our 
             ! heights calculated will be zero. In the IF-statement we could use .LT.
             ! min_stomate but that would mean that at right after coppicing the
             ! canopy will be 1 m heigh but at as soon as we grow a bit of sapwood
             ! the height could shrink again because we have too much sapwood to satisfy
             ! the IF-statement but not enough to make a canopy that is higher than 1 m. 
             ! For that reason we used .LT. un. In that way we will not have shrinking
             ! canopies but we might prescribe the canopy height in this calculation
             ! for a much longer period.
             IF (SUM(circ_class_biomass(ipts,ivm,:,isapabove,icarbon)) < min_stomate .AND. &
                 SUM(circ_class_biomass(ipts,ivm,:,isapbelow,icarbon)) > min_stomate) THEN

                !CALL ipslerr_p(3,'stomate_laieff.f90', 'The solution in the code has not yet been tested', &
                !     'Remove this stop and test the solution given in the code', '')
                
                ! Prescribe a small canopy. 
                max_height = un
                min_height = 0.0001_r_std

             END IF

             ! Debug
             IF(printlev_loc.GT.4)THEN
               IF((ivm == test_pft) .AND. (ipts == test_grid))THEN
                   WRITE(numout,*) 'This is the largest height of trees in z_level_photo'
                   WRITE(numout,*) 'max_height, min_height: ',max_height,min_height
                   WRITE(numout,*) 'heights: ',heights(:)
                   WRITE(numout,*) 'cn_dia: ',cn_dia(:)
                   DO il=1,nlevels
                      WRITE(numout,*) 'il, ivm: ',il,ivm
                      WRITE(numout,*) 'min_level_sep(ivm): ',&
                        min_level_sep(ivm)
                      IF(il .NE. nlevels) &
                           WRITE(numout,*) 'z_level(il+1),z_level(il):',z_level(il+1),z_level(il)
                   ENDDO
                ENDIF
             ENDIF
             !-

          ELSE

             ! Height for grasses and croplands in line with stomate_growth_fun_all.f90
             max_height = biomass_to_lai(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)* &
                circ_class_n(ipts,ivm,:)),ivm) * lai_to_height(ivm)

             ! arbitrary. Don't want it equal to zero because
             ! it won't enter the correct loop below.
             min_height = 0.0001_r_std 

             ! Debug
             IF(printlev_loc.GT.4)THEN
                IF((ivm == test_pft) .AND. (ipts == test_grid))THEN
                   WRITE(numout,*) 'This is the largest height of grasses/crops in z_level_photo'
                   WRITE(numout,*) 'max_height, min_height: ',max_height,min_height
                   WRITE(numout,*) 'lai_to_height: ',lai_to_height(ivm)
                ENDIF
             ENDIF
             !-

          ENDIF
  
          ! Loop over the levels used by the energy budget.  For each level, there are
          ! six possibilities, based on the location of the canopy vis-a-vis the location
          ! of this level height and the level above (if there is no level above, we take
          ! a boundary that is equal to a very big number).
          ! The word "slice" refers to the space between this level and the level above.
          ! 1) Only the bottom part of the canopy is in this slice.
          ! 2) Only the top part of the canopy is in this slice.
          ! 3) The whole canopy is below this slice.
          ! 4) The whole canopy is above this slice.
          ! 5) The whole canopy is inside this slice.
          ! 6) The top of the canopy is above the slice, the bottom is below the
          !    slice, and we just have the middle part in the slice.
          ! We adjust the photosynthesis levels in this slice to cover as much as the canopy as possible, assuming
          ! a mininum distince between the levels.  This is done to avoid having a lot of levels with
          ! very little LAI, which could cause a problem in the photosyntheis routines.
          ! If there is too much LAI, the absorbed light saturates and we lose
          ! sensitivity in the photosynthesis routines.  0.5 seems to be a good LAI
          ! to have in a level, but this is just a feeling.
          DO il=1,nlevels

             lower_level_boundary=z_level(il)

             ! The lowest photosynthesis level should always be equal to
             ! the lower_level_boundary.  In this way we make sure that
             ! we account for all the vegetation in this level.  We need
             ! to make sure we don't overwrite this below.  This means that
             ! the highest PS level in this level should always be at least
             ! min_level_sep from the upper_boundary, and that instead of
             ! nlai levels to play with, we only have nlai-1.  This
             ! should not cause any divide by zero errors since there is
             ! a loop at the beginning of the subroutine to handle the trivial
             ! case where nlai=1.
             z_level_photo(ipts,ivm,nlai*(il-1)+1)=lower_level_boundary

             IF(il .NE. nlevels)THEN
                upper_level_boundary=z_level(il+1)
             ELSE
                ! This number is arbitary, but it just needs to be larger than
                ! the largest possible vegetation height.  10000 meters is a very
                ! tall tree.
                upper_level_boundary=10000.0_r_std
             ENDIF

             ! Visit all the possible locations of the canopy with respect to
             ! the level boundaries.
             IF((min_height .GT. lower_level_boundary .AND. &
                  min_height .LT. upper_level_boundary ) .AND. &
                  (max_height .GT. upper_level_boundary))THEN 
                ! 1) This slice only has the bottom part of the canopy in it.
                ! Ideally, we would want the second lowest PS level to be just
                ! above the bottom of the canopy, and the highest PS level
                ! to be just below the upper boundary.  If there is not
                ! enough of the canopy to do this and maintain our
                ! minimum level spacing, we just space the levels
                ! with this minimum starting from just below the upper boundary.
                ideal_level_thickness = ((upper_level_boundary-min_level_sep(ivm))-&
                     (min_height+min_level_sep(ivm)))&
                     /REAL(nlai-2,r_std)

                IF(ideal_level_thickness .GT. min_level_sep(ivm))THEN
                   ! It worked, the levels are thick enough.
                   DO il_photo=nlai,2,-1
                      z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                           (upper_level_boundary-min_level_sep(ivm))-&
                           ideal_level_thickness*REAL(nlai-il_photo,r_std)
                   ENDDO

                ELSE

                   ! Use the minimum level thickness starting from the 
                   ! upper boundary.
                   DO il_photo=nlai,2,-1
                      z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                           (upper_level_boundary-min_level_sep(ivm))-&
                           min_level_sep(ivm)*REAL(nlai-il_photo,r_std)
                   ENDDO

                ENDIF

             ELSEIF(min_height .LT. lower_level_boundary .AND. &
                  ((max_height .GT. lower_level_boundary) .AND. &
                  (max_height .LT. upper_level_boundary)))THEN
                ! 2) Only the top half of the canopy is in this slice.  We need
                ! to do the same thing as case (1), but in the other direction.
                ideal_level_thickness = ((max_height-min_level_sep(ivm))-&
                     (lower_level_boundary+min_level_sep(ivm)))&
                     /REAL(nlai-2,r_std)

                IF(ideal_level_thickness .GT. min_level_sep(ivm))THEN
                   ! It worked, the levels are thick enough.
                   DO il_photo=2,nlai
                      z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                           (lower_level_boundary+min_level_sep(ivm))+&
                           ideal_level_thickness*REAL(il_photo-2,r_std)
                   ENDDO

                ELSE

                   ! Use the minimum level thickness starting from the 
                   ! upper boundary.
                   DO il_photo=2,nlai
                      z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                           (lower_level_boundary+min_level_sep(ivm))+&
                           min_level_sep(ivm)*REAL(il_photo-2,r_std)
                   ENDDO

                ENDIF

             ELSEIF((min_height .GT. upper_level_boundary) .OR. &
                  (max_height .LT. lower_level_boundary) )THEN
                ! (3) and (4) The whole canopy is below the slice, or the whole canopy
                ! is above this slice.  Whatever we do does not
                ! matter, since there is no LAI in the slice and therefore
                ! there will be no PS.  
                DO il_photo=2,nlai
                   z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=lower_level_boundary+&
                        min_level_sep(ivm)*REAL(il_photo-1,r_std)
                ENDDO

             ELSEIF((min_height .GT. lower_level_boundary) .AND. &
                  (max_height .LT. upper_level_boundary) )THEN
                ! (5) The whole canopy is inside this slice.  This is the trickiest
                ! case, and unfortunately the most common.  In fact, it's the only
                ! case which can exist with a single level energy budget.
                ideal_level_thickness = ((max_height-min_level_sep(ivm))-&
                     (min_height+min_level_sep(ivm)))&
                     /REAL(nlai-2,r_std)

                IF(ideal_level_thickness .GT. min_level_sep(ivm))THEN
                   ! It worked, the levels are thick enough.
                   DO il_photo=2,nlai
                      z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                           (min_height+min_level_sep(ivm))+&
                           ideal_level_thickness*REAL(il_photo-2,r_std)
                   ENDDO

                ELSE

                   ! Here we have a block of min_level_sep levels.  Let's place the 
                   ! center of this block in the center of our canopy, and check to see
                   ! if the ends poke past the boundaries of this slice.
                   block_thickness = REAL(nlai-2,r_std) * min_level_sep(ivm)
                   canopy_thickness = (max_height - min_height )
                   canopy_center = min_height + canopy_thickness / deux

                   IF( ((canopy_center + block_thickness/deux) .LT. &
                        (upper_level_boundary - min_level_sep(ivm))) .AND. &
                        (canopy_center - block_thickness/deux) .GT. &
                        (lower_level_boundary + min_level_sep(ivm)))THEN
                      ! The block center can be placed at the center of the canopy, which means that
                      ! the bottom of the block starts at canopy_center - block_thickness/2

                      DO il_photo=2,nlai
                         z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                              canopy_center - block_thickness/deux+&
                              min_level_sep(ivm)*REAL(il_photo-2,r_std)
                      ENDDO

                   ELSEIF( (canopy_center - block_thickness/deux) .LE. &
                        (lower_level_boundary + min_level_sep(ivm)))THEN
                      ! we have enough space up, but not down, so move the whole block upwards.
                      ! All we have to do is start the block above the lower boundary level.

                      DO il_photo=2,nlai
                         z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                              lower_level_boundary + &
                              min_level_sep(ivm)*REAL(il_photo-1,r_std)
                      ENDDO

                   ELSEIF( (canopy_center + block_thickness/deux) .GE. &
                        (upper_level_boundary - min_level_sep(ivm)))THEN
                      ! we have enough space down, but not up, so move the whole block downwards
                      ! All we need to do is start the block just below the upper boundary.

                      DO il_photo=nlai,2,-1
                         z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                              upper_level_boundary - &
                              min_level_sep(ivm)*REAL(nlai-il_photo+1,r_std)
                      ENDDO

                   ELSE

                      ! We should not be here!
                      WRITE(numout,*) 'ERROR: In calculate_z_photo_level, we have reached a '
                      WRITE(numout,*) 'ERROR: case that was not foreseen!'
                      WRITE(numout,*) 'ipts,ivm,il: ',ipts,ivm,il
                      WRITE(numout,*) 'canopy_center: ',canopy_center
                      WRITE(numout,*) 'block_thickness: ',block_thickness
                      WRITE(numout,*) 'upper_level_boundary: ',upper_level_boundary
                      WRITE(numout,*) 'lower_level_boundary: ',lower_level_boundary
                      WRITE(numout,*) 'nlai: ',nlai
                      CALL ipslerr_p (3,'calculate_z_level_photo', 'Unforeseen case.','','')

                   ENDIF

                ENDIF

             ELSEIF((min_height .LT. lower_level_boundary) .AND. &
                  (max_height .GT. upper_level_boundary) )THEN

                ! (6) This slice is completely filled by the canopy.  This is easy,
                ! we just evenly space the levels in this slice.
                ideal_level_thickness = ((upper_level_boundary-min_level_sep(ivm))-&
                     (lower_level_boundary+min_level_sep(ivm)))&
                     /REAL(nlai-2,r_std)

                DO il_photo=2,nlai
                   z_level_photo(ipts,ivm,nlai*(il-1)+il_photo)=&
                        (lower_level_boundary+min_level_sep(ivm))+&
                        ideal_level_thickness*REAL(il_photo-2,r_std)
                ENDDO
          
             ELSE

                ! Should not be here!  This is a case I haven't thought of.
                WRITE(numout,*) 'Should not be here in find photosynthesis levels!'
                WRITE(numout,*) 'ipts,ivm,il: ',ipts,ivm,il
                WRITE(numout,*) 'max_height: ',max_height
                WRITE(numout,*) 'min_height: ',min_height
                WRITE(numout,*) 'upper_level_boundary: ',upper_level_boundary
                WRITE(numout,*) 'lower_level_boundary: ',lower_level_boundary
                CALL ipslerr_p (3,'calculate_z_level_photo', 'Unforeseen case in ps levels.','','')
                
             ENDIF 

          ENDDO

          ! Debug
          IF(printlev_loc >= 4 .AND. ivm == test_pft .AND. ipts == test_grid)THEN
             WRITE(numout,*) 'Printing the level heights used by photosynthesis.'
             WRITE(numout,*) 'ipts,ivm,max_height: ',ipts,ivm,max_height
             DO il=1,nlevels_tot
                WRITE(numout,*) 'il,z_level_photo(ipts,ivm,il): ',&
                     il,z_level_photo(ipts,ivm,il)
             ENDDO
          ENDIF
          !-

       ENDDO
    ENDDO

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving calculate_z_level_photo'

  END SUBROUTINE calculate_z_level_photo

!! ================================================================================================================================
!! SUBROUTINE 	: combine_lai_levels
!!
!>\BRIEF        Collapses the LAI per photosynthesis level into an array
!!              which only has the LAI per energy budget level.
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai_per_level_temp
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE combine_lai_levels(lai_per_level, lai_per_level_temp)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    REAL(r_std),DIMENSION(:,:,:),INTENT(IN)        :: lai_per_level        !! The LAI per photosynthesis level
                                                                           !! @tex $(m^{2} m^{-2})$


    !! 0.2 Output variables
    REAL(r_std),DIMENSION(:,:,:),INTENT(OUT)       :: lai_per_level_temp   !! The LAI per energy budget level
                                                                           !! @tex $(m^{2} m^{-2})$

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                        :: il, ilevel           !! indices                        
!_ ================================================================================================================================
    IF (firstcall_laieff) CALL ipslerr_p (3,'combine_lai_levels', 'Initialization is not done before first call','','')

    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering combine_lai_levels'

    lai_per_level_temp(:,:,:)=zero

    DO ilevel=1,nlevels_tot
       il=CEILING(REAL(ilevel)/REAL(nlai))
       lai_per_level_temp(:,:,il)=lai_per_level_temp(:,:,il)+lai_per_level(:,:,ilevel)
    ENDDO

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving combine_lai_levels'

  END SUBROUTINE combine_lai_levels

!! ================================================================================================================================
!! SUBROUTINE 	: find_lai_per_level
!!
!>\BRIEF         Calculates the LAI per vertical level of the canopy.
!!
!! DESCRIPTION  :  Given a total biomass, number of individuals, size of individuals,
!!                 and set of vertical layer boundaries, this routine calculates the
!!                 the leaf area index (LAI) found between the layer boundaries.  It
!!                 first makes a call to the stand structure in order to find out some
!!                 canopy properties like height and canopy diameter, and then uses geometry
!!                 to determine the amount of LAI present.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai_per_level
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE find_lai_per_level(npts, z_level_photo, &
       circ_class_biomass, circ_class_n, lai_per_level, max_height_store)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(IN)                         :: npts               !! Domain size - number of pixels (unitless)
    REAL(r_std),DIMENSION(:,:,:),INTENT(IN)           :: z_level_photo      !! The height of the bottom of each of our levels @tex $(m)$ @endtex
                                                                            !! Note that level 1 is closest to the ground
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN)     :: circ_class_biomass !! Biomass of the componets of the model  
                                                                            !! tree within a circumference
                                                                            !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)         :: circ_class_n       !! Number of trees within each circ class
                                                                            !! class @tex $(m^{-2})$ @endtex
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(OUT)        :: lai_per_level      !! This is the LAI per vertical level
                                                                            !! @tex $(m^{2} m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(OUT)          :: max_height_store


    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER                                            :: ipts, ivm, il, il_photo, &
                                                          icir                  !! indices
    REAL(r_std)                                        :: max_height            !! The height of the tallest tree
                                                                                !! @tex $(m)$ @endtex
    REAL(r_std)                                        :: ind_loc               !! The number of individuals
                                                                                !! @tex $(m^{-2})$ @endtex
    REAL(r_std)                                        :: upper_boundary        !! The top of the current level 
                                                                                !! @tex $(m)$ @endtex
    REAL(r_std)                                        :: lower_boundary        !! The bottom of the current level 
                                                                                !! @tex $(m)$ @endtex
    REAL(r_std)                                        :: lai_temp              !! LAI of the current PFT/grid square
                                                                                !! @tex $(m^{2} m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,ndist_types)  :: values                !! An array which holds various canopy parameters
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot)        :: PartialCrownVolumeL   !! The fraction of the crown volume found
                                                                                !! in this level averaged over all circ classes
                                                                                !! @tex $(m^{3})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_tot)  :: PartialCrownVolume_h  !! The fraction of the crown volume found
                                                                                !! in this level for each circ class
                                                                                !! @tex $(m^{3})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm)                    :: CrownVolumeL          !! The average crown volume of the stand.
                                                                                !! @tex $(m^3)$ @endtex
    REAL(r_std),DIMENSION(npts,nvm)                    :: favd                  !! Foliage area volume density
                                                                                !! @tex $(m^2 (leaf) m^{-3})$ @endtex

!_ ================================================================================================================================

    IF (firstcall_laieff) CALL ipslerr_p (3,'find_lai_per_level', 'Initialization is not done before first call','','')

    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering find_lai_per_level' 

    CALL derive_biomass_quantities(npts, nvm, ncirc, circ_class_n, &
         circ_class_biomass, values)

    PartialCrownVolumeL(:,:,:)=zero
    PartialCrownVolume_h(:,:,:,:)=zero
    CrownVolumeL(:,:)=zero
    favd(:,:) = zero    

    DO ipts=1,npts
       DO ivm=1,nvm

          ! Skip all of this if we have no leaf biomass
          IF(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)) .LE. min_stomate)THEN
             lai_per_level(ipts,ivm,:)=zero
             CYCLE
          ENDIF

          ind_loc=SUM(circ_class_n(ipts,ivm,:))

          IF(ind_loc .LT. min_stomate)THEN
             WRITE(numout,*) 'How do we have biomass but no individuals?'
             WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
             WRITE(numout,*) 'biomass(ipts,ivm,ileaf,icarbon),ind_loc: ',&
                  SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)),ind_loc
             CALL ipslerr_p (3,'find_lai_per_level', 'Biomass but no individuals.','','')
          ENDIF

          CALL calculate_favd(ipts,ivm, values(ipts,ivm,:,icnvol), circ_class_biomass, &
               circ_class_n, favd)

          ! Find the LAI of this PFT.
          lai_temp = biomass_to_lai(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)*circ_class_n(ipts,ivm,:)), ivm)
        

          ! This is done differently for forests and grass/crops, since grass
          ! and crops will not have a canopy volume.  We just assume a solid
          ! block of vegetation with a given height.

          IF(is_tree(ivm))THEN

             ! Debug
             IF(printlev_loc.GT.4 .AND. test_pft == ivm .AND. ipts == test_grid)&
                  WRITE(6,*) 'On tree ',ivm,lai_temp,favd(ipts,ivm)
             !-

             max_height = MAXVAL(values(ipts,ivm,:,iheight))

             ! Total crown volume weighted by the number of trees in each circ class
             ! We need to integrate over the height distribution. This is just a
             ! simple rectangular integration, and it gives us the average crown
             ! volume for the whole stand.
             DO icir=1,ncirc
                CrownVolumeL(ipts,ivm)=CrownVolumeL(ipts,ivm)+&
                     values(ipts,ivm,icir,icnvol)*&
                     circ_class_n(ipts,ivm,icir)/ind_loc
                IF(test_pft == ivm .AND. printlev_loc.GT.4 .AND. ipts == test_grid)THEN
                   WRITE(numout,*) 'Volume for icir: ',icir,values(ipts,ivm,icir,icnvol)
                ENDIF
             ENDDO

             ! Now I need to find the average crown volume per layer.  This is not the
             ! same thing as taking the average crown volume and dividing it up
             ! over the layers, since the crowns are at different heights.
             DO il=1,nlevels_tot
                
                IF(il .NE. nlevels_tot)THEN
                   upper_boundary=z_level_photo(ipts,ivm,il+1)
                ELSE
                   ! This height is arbitrary.  It just has to include the whole canopy.
                   upper_boundary=z_level_photo(ipts,ivm,nlevels_tot)+max_height
                ENDIF
                lower_boundary=z_level_photo(ipts,ivm,il)

                ! For this level, how much of each canopy is here?  We will also
                ! integrate this value over the whole height distribution, thus
                ! finding the average canopy volume in this layer.
                DO icir=1,ncirc

                   PartialCrownVolume_h(ipts,ivm,icir,il)=&
                        partial_spheroid_vol(upper_boundary, lower_boundary, &
                        values(ipts,ivm,icir,icndiaver)/deux, &
                        values(ipts,ivm,icir,icndiahor)/deux, &
                        values(ipts,ivm,icir,iheight))

                   ! Debug
                   IF(printlev_loc.GT.4 .AND. test_pft == ivm .AND. test_grid == ipts)THEN
                      WRITE(numout,*) 'ipts,ivm,icir,il: ',ipts,ivm,icir,il
                      WRITE(numout,*) 'PartialCrownVolume_h(ipts,ivm,icir,il): ',&
                           PartialCrownVolume_h(ipts,ivm,icir,il)
                      WRITE(numout,*) 'crown diameter: ', values(ipts,ivm,icir,icndiaver)
                      WRITE(numout,*) 'tree height: ', values(ipts,ivm,icir,iheight)
                      WRITE(numout,*) 'circ_class_n(ipts,ivm,icir)',circ_class_n(ipts,ivm,icir)
                   ENDIF
                   !-

                    ! It's important to note here that the units of PartialCrownVolume_h
                    ! are m**3/tree, since the canopy properties that we passed to
                    ! partial_spheroid_vol are also per tree.  We need to get it in
                    ! m**3/m**2 (land area).
                    PartialCrownVolumeL(ipts,ivm,il)=PartialCrownVolumeL(ipts,ivm,il)+&
                         PartialCrownVolume_h(ipts,ivm,icir,il)*&
                         circ_class_n(ipts,ivm,icir)

                END DO

                ! Now I have the average crown volume per level.  I also have
                ! the foliage area volume density.  The amount of LAI in this
                ! level is just the combination of the two.
                lai_per_level(ipts,ivm,il)=PartialCrownVolumeL(ipts,ivm,il)*&
                     favd(ipts,ivm)

             ENDDO

              IF(printlev_loc.GT.4 .AND. ipts == test_grid .AND. ivm == test_pft)THEN
                 WRITE(numout,*) 'CrownVolumeL: ',CrownVolumeL(ipts,ivm)
              ENDIF

          ELSE

             ! For grasses and crops.

             ! This is not ideal, since we already had to calculate lai_temp 
             ! and max_height to find the favd.  If there is a speed problem,
             ! we can change it.

             ! What is the height of our grass/crop?
       
             ! This converts LAI to a height.  I took this from functional_allocation.
             max_height = lai_temp * lai_to_height(ivm)


             ! Find out how much of this vegetation cube is in this layer
             ! These calculations here include some implict factors of 
             ! m**2/m**2, so be aware of that when calculating the levels.
             ! All heights are multiplied by 1x1 m to find a volume, but
             ! that is not explictly shown.
             DO il=1,nlevels_tot

                IF(il .NE. nlevels_tot)THEN
                   upper_boundary=z_level_photo(ipts,ivm,il+1)
                ELSE
                   ! This height is arbitrary.  It just has to include the whole canopy.
                   upper_boundary=z_level_photo(ipts,ivm,nlevels_tot)+max_height
                ENDIF
                lower_boundary=z_level_photo(ipts,ivm,il)

                IF(max_height .GE. upper_boundary)THEN
                   ! This whole level is full of vegetation.
                   lai_per_level(ipts,ivm,il)=&
                        (upper_boundary-lower_boundary)*favd(ipts,ivm)

                ELSEIF(max_height .LE. lower_boundary)THEN

                   ! There is no vegetation in this level.
                   lai_per_level(ipts,ivm,il)=zero

                ELSE
                   
                   ! This level is partially full.
                   lai_per_level(ipts,ivm,il)=&
                        (max_height-lower_boundary)*favd(ipts,ivm)

                ENDIF ! checking where max_height is

             ENDDO

          ENDIF

          IF(test_pft == ivm .AND. test_grid == ipts .AND. printlev_loc.GT.4)THEN
             WRITE(numout,*) 'LAI per level: ',ipts,ivm
             DO il=nlevels_tot,1,-1
                WRITE(numout,*) 'il,lai_per_level(ipts,ivm,il): ',&
                     il,lai_per_level(ipts,ivm,il)
             ENDDO
          ENDIF

          ! Test to see if our LAI is equal to the sum over all the levels.
          IF( ABS(lai_temp - SUM(lai_per_level(ipts,ivm,:))) &
               .GT. 1e-5)THEN
             WRITE(numout,*) '!*****************************************************'
             WRITE(numout,*) 'Problem in LAI levels! '
             WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
             WRITE(numout,*) 'lai_temp, SUM(lai_per_level(ipts,ivm,:)): ',&
                  lai_temp,SUM(lai_per_level(ipts,ivm,:))
             DO il=nlevels_tot,1,-1
                WRITE(numout,*) 'il,lai_per_level(ipts,ivm,il): ',&
                     il,lai_per_level(ipts,ivm,il)
             ENDDO

             ! There can be issues using the functions in the print statements.  If a print 
             ! statement in a function is activated when the function call itself is in a print
             ! statement, there is a crash.
             WRITE(numout,*) 'Check whether canopy dimensions exceed the height '
             WRITE(numout,*) 'diameter, ', values(ipts,ivm,:,idiameter)
             WRITE(numout,*) 'height: ', values(ipts,ivm,:,iheight)
             WRITE(numout,*) 'crown area: ', values(ipts,ivm,:,icnarea)
             WRITE(numout,*) 'crown diameter - ver: ', values(ipts,ivm,:,icndiaver)
             WRITE(numout,*) 'crown diameter - hor: ', values(ipts,ivm,:,icndiahor)
             WRITE(numout,*) 'crown volume: ', values(ipts,ivm,:,icnvol)
             WRITE(numout,*) 'crown diameter used in this routine: ', (6./pi*values(ipts,ivm,:,icnvol))**(1./3.) 
             WRITE(numout,*) '!*****************************************************'
             CALL ipslerr_p (3,'find_lai_per_level', 'Inconsistency with LAI in levels.','','')
          ENDIF ! check to see if the LAI per level is consistent with the total LAI

          max_height_store(ipts, ivm) = max_height


       ENDDO ! loop over PFTs
    ENDDO ! loop over points

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving find_lai_per_level'

  END SUBROUTINE find_lai_per_level


!! ================================================================================================================================
!! SUBROUTINE 	: fitting_laieff
!!
!>\BRIEF          Fits the effective LAI as a function of the solar angle    
!!
!! DESCRIPTION  :  Calculation of the effective LAI is a very expensive part of the
!!                 model, due to the loops over points, PFTs, circ classes, and levels,
!!                 in addition to using trig functions.  We need this value at
!!                 every time step (up to 48 times a day).  Thankfully, the function
!!                 doesn't change a lot, in particular when the sun is high in the sky.
!!                 Therefore, we will calculate a few points exactly, and then fit
!!                 a function to these points to predict the value of the LAIeff at any
!!                 solar angle.  This method seems to work well for low solar zenith
!!                 angles, which are the most important since they will have the most
!!                 incoming solar radiation.
!!
!!    NOTE:
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::laieff_fit
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE fitting_laieff(npts, z_array_out, circ_class_biomass, &
       circ_class_n, veget_max, lai_per_level, laieff_fit)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(IN)                              :: npts                !! Domain size - number of pixels (unitless)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(IN)              :: z_array_out         !! the height of the bottom of each of our levels 
                                                                                  !! @tex $(m)$ @endtex
                                                                                  !! note that level 1 is closest to the ground
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN)          :: circ_class_biomass  !! Biomass of the componets of the model  
                                                                                  !! tree within a circumference class
                                                                                  !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)              :: circ_class_n        !! Number of trees within each circ class
                                                                                  !! @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(IN)                  :: veget_max           !! Maximum fraction of vegetation type including 
                                                                                  !! non-biological fraction (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)              :: lai_per_level       !! This is the LAI per vertical level
                                                                                  !! @tex $(m^{2} m^{-2})$
 

    !! 0.2 Output variables
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(OUT)        :: laieff_fit          !! Fitted parameters for the effective LAI


    !! 0.3 Modified variables

    !! 0.4 Local variables

    !+++CHECK+++
    ! Idealy, the next two values would be externalised.
    INTEGER, PARAMETER                                     :: nangles=8           !! The number of angles to fit the parameterization
                                                                                  !! of the effective LAI. 

    REAL(r_std), DIMENSION(nangles)                        :: param_angles        !! the angles to
    !+++++++++++

    REAL(r_std), DIMENSION(nlevels_tot,npts,nvm,nangles)   :: laieff              !! Leaf Area Index Effective converts 3D lai into 
                                                                                  !! 1D lai for two stream radiation transfer model
                                                                                  !! @tex $(m^{2} m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                           :: solar_angle         !! The solar zenith angle of every grid square 
                                                                                  !! @tex $(rad)$ @endtex 
    INTEGER(i_std)                                         :: ilevel, ipts, icir  !! index (unitless)
    INTEGER(i_std)                                         :: ivm, iangle, jangle !! index (unitless)

    INTEGER(i_std)                                         :: nzero               !! Counter to keep track of the number of
                                                                                  !! LAIeff with value zero
    INTEGER(i_std)                                         :: temp_index          !! An index which allows us to fit only non-zero
                                                                                  !! values
    REAL(r_std), DIMENSION(npts)                           :: test_angles         !! Converts the param_angles (degrees) into radians
                                                                                  !! @tex $(rad)$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot,nangles)    :: Pgap_cumul          !! The probability of finding a gap in the
                                                                                  !! in canopy from the top of the canopy
                                                                                  !! to a given level.
                                                                                  !! (unitless, between 0-1)
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot,nangles)    :: Pgap_perlevel       !! The probability of finding a gap in the
                                                                                  !! in the canopy between a level
                                                                                  !! and the level immediately above.
                                                                                  !! (unitless, between 0-1)
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot)            :: temp                !! temporary variable for history files
    REAL(r_std), DIMENSION(npts)                           :: fedge_pgap          !! Edge effect weight on canopy porosity (0-1)
                                                                                  !! See MODULE_DESIGN_AED_LIGHT.md §4.2

!_ ================================================================================================================================

    IF (firstcall_laieff) CALL ipslerr_p (3,'fitting_laieff', 'Initialization is not done before first call','','')

    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering fitting_laieff'

    ! Debug
    IF (printlev_loc.GT.4) THEN
       WRITE(numout, *) '130416 nlevels_tot is: ', nlevels_tot
       WRITE(numout, *) '130416 nlai is: ', nlai
    ENDIF
    !-

    !! 1. Initialize check for surface area conservation
    !  Veget_max is a INTENT(in) variable and can therefore
    !  not be changed during the course of this subroutine
    !  No need to check whether the subroutine preserves the
    !  total surface area of the pixel.

    !! 2. Calculations
    ! We want to parameterize the effective LAI as a function of
    ! the solar zenith angle for every level, grid square, and
    ! PFT type.  Let's choose a series of angles to use.
    ! There are somtimes numerical issues with using exactly
    ! 90 degrees or exactly 0 degrees.

!    param_angles = (/ 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0 /)
    ! Using high SZA (near the horizon) seems to cause problems
    ! with the fitting, as strange behavior happens there.  Since there
    ! is not a lot of light, let's use smaller values.
    param_angles = (/ 5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0 /)

    ! Right now, I'm doing this the brute force way.  It will probably
    ! need to be refined later.
    DO iangle=1,nangles

       test_angles(:)=COS(param_angles(iangle)/180.0*pi)
       
       ! Debug
       IF(printlev_loc.GT.4)&
            WRITE(numout,*) 'Calculating for test angle: ',param_angles(iangle)
       !-

       CALL effective_lai(npts, nlevels_tot, z_array_out, circ_class_biomass, &
            circ_class_n, veget_max, test_angles, lai_per_level, &
            laieff(:,:,:,iangle), Pgap_cumul_out=Pgap_cumul(:,:,:,iangle), &
            Pgap_perlevel_out=Pgap_perlevel(:,:,:,iangle))

       ! Debug
       IF(printlev_loc.GT.4)THEN
          DO ipts=1,npts
             DO ivm=1,nvm
                IF(ipts == test_grid .AND. ivm == test_pft)THEN
                   WRITE(numout,*) 'LAIEFF for angle ',param_angles(iangle)
                   DO ilevel=nlevels_tot,1,-1
                      WRITE(numout,*) ilevel,laieff(ilevel,ipts,ivm,iangle)
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !-

    ENDDO

    !! 2bis. Non-directional canopy edge light effect (Marie 2026 / MODULE_DESIGN_AED_LIGHT.md §4.2)
    ! Apply post-effective_lai modulation of Pgap_cumul / Pgap_perlevel / laieff
    ! by fedge_pgap(AED_light). Cascades downstream through laieff_fit -> Diffuco -> GPP.
    ! The companion site that drives the albedo cascade is slowproc_veget (slowproc.f90:3770).
    ! AED_light comes from stomate_data.calculate_aed (called per-day from stomate.f90).
    IF (ok_aed_light .AND. ALLOCATED(AED_light)) THEN
       fedge_pgap(:) = 1.0_r_std - ((AED_light(:) - edge_distance_light/2._r_std)**2) &
                                  / (AED_light(:)**2)

       DO iangle = 1, nangles
          DO ilevel = 1, nlevels_tot
             DO ipts = 1, npts
                Pgap_cumul(ipts,:,ilevel,iangle) =                              &
                   Pgap_cumul(ipts,:,ilevel,iangle) *                           &
                   (1.0_r_std + pgap_edge_factor * fedge_pgap(ipts))
                Pgap_perlevel(ipts,:,ilevel,iangle) =                           &
                   Pgap_perlevel(ipts,:,ilevel,iangle) *                        &
                   (1.0_r_std + pgap_edge_factor * fedge_pgap(ipts))
                laieff(ilevel,ipts,:,iangle) =                                  &
                   laieff(ilevel,ipts,:,iangle) *                               &
                   (1.0_r_std + laieff_edge_factor * fedge_pgap(ipts))
             ENDDO
          ENDDO
       ENDDO

       Pgap_cumul(:,:,:,:)    = MIN(1.0_r_std, MAX(0.0_r_std, Pgap_cumul(:,:,:,:)))
       Pgap_perlevel(:,:,:,:) = MIN(1.0_r_std, MAX(0.0_r_std, Pgap_perlevel(:,:,:,:)))
       laieff(:,:,:,:)        = MAX(0.0_r_std, laieff(:,:,:,:))

       CALL xios_orchidee_send_field("AED_LIGHT",  AED_light)
       CALL xios_orchidee_send_field("FEDGE_PGAP", fedge_pgap)
    ENDIF

    ! Write the gap probability profile throughout the canopy for a solar
    ! zenith angle of 5 degrees
    CALL xios_orchidee_send_field("P_GAP_LEVEL",Pgap_perlevel(:,:,:,1))
    CALL xios_orchidee_send_field("P_GAP_CUM",Pgap_cumul(:,:,:,1))

    param_angles(:)=COS(param_angles(:)/180.0_r_std*pi)

    ! Now we have all the data points, so we can fit the function.
    DO ipts=1,npts
       DO ivm=1,nvm
          DO ilevel=1,nlevels_tot

             ! The 0.001 is completely arbitrary here.
             IF(SUM(laieff(ilevel,ipts,ivm,:)) .GT. 0.001_r_std )THEN

                ! We have a choice to make.  Sometimes the laieff for highest zenith angles
                ! will be zero, even if the previous points were not.  This seems to happen
                ! when no light gets through, which is normal for lower vegetation
                ! levels as the sun gets closer to the horizon.  The pattern of points is
                ! something like 1,1,1,1,1,1.9,4.0,0.0.  This causes a problem for
                ! fitting, and we don't care too much about this largest angle.  If we
                ! can fit the first points well, and the last point is greater than 4.0, 
                ! we will probably be okay since there is not a lot of sunlight wen the
                ! sun is near the horizon in any case.  So let's run a test.
                IF(laieff(ilevel,ipts,ivm,nangles) .LT. min_stomate)THEN

                   ! We know that the sum of all elements is greater than our cutoff
                   ! above, so we need to find the final substantial element.
                   ! Count backwards to see how many zero elements we have
                   ! at the end.
                   nzero=1
                   DO iangle=nangles-1,1,-1
                      IF(laieff(ilevel,ipts,ivm,iangle) .LT. min_stomate)THEN
                         nzero=nzero+1
                      ELSE
                         EXIT
                      ENDIF
                   ENDDO
                   
                   IF(nzero .GE. nangles)THEN
                      WRITE(numout,*) 'ERROR: Something went wrong in fitting the laieff!'
                      WRITE(numout,*) 'ipts,ivm,ilevel: ',ipts,ivm,ilevel
                      WRITE(numout,*) 'nzero,nangles: ',nzero,nangles
                      WRITE(numout,*) 'laieff(ilevel,ipts,ivm,:): ',laieff(ilevel,ipts,ivm,:)
                      CALL ipslerr_p (3,'fitting_laieff', 'All zero values in the fit.','','')
                   ENDIF

                   ! This three is also arbitrary.
                   IF(nzero .GE. 3)THEN
                      WRITE(numout,*) 'ERROR: A surprising number of zero values in laieff!'
                      WRITE(numout,*) 'ipts,ivm,ilevel: ',ipts,ivm,ilevel
                      WRITE(numout,*) 'nzero,nangles: ',nzero,nangles
                      WRITE(numout,*) 'laieff(ilevel,ipts,ivm,:): ',laieff(ilevel,ipts,ivm,:)
                      CALL ipslerr_p (3,'fitting_laieff', 'Too many zero values in the fit.','','')
                   ENDIF

                   ! Only use the first non-zero values for the fitting
                   temp_index=nangles-nzero
                   CALL fit_laieff_to_points(temp_index,param_angles(1:temp_index),&
                        laieff(ilevel,ipts,ivm,1:temp_index), &
                        laieff_fit(ipts,ivm,ilevel),ipts,ivm,ilevel)
                ELSE
                   CALL fit_laieff_to_points(nangles,param_angles,laieff(ilevel,ipts,ivm,:), &
                        laieff_fit(ipts,ivm,ilevel),ipts,ivm,ilevel)
                ENDIF

             ELSE

                ! Not enough vegetation to fit it.  Everything is zero.
                laieff_fit(ipts,ivm,ilevel)%a=zero
                laieff_fit(ipts,ivm,ilevel)%b=zero
                laieff_fit(ipts,ivm,ilevel)%c=zero
                laieff_fit(ipts,ivm,ilevel)%d=zero
                laieff_fit(ipts,ivm,ilevel)%e=zero
             ENDIF
          END DO
       ENDDO
    ENDDO

    ! Now we have all the data points, so we can fit the function.
    DO ipts=1,npts
       DO ivm=1,nvm
          DO ilevel=1,nlevels_tot
                  IF(laieff_fit(ipts,ivm,ilevel)%a .LT. 0.0d0)THEN
                     ! temporary, during merge, to try and eliminate negative values
                     IF(printlev_loc.GT.4) &
                          WRITE(numout, *) '130416c laieff is at this point less than zero' 
                     laieff_fit(ipts,ivm,ilevel)%a = 0.0d0
                  END IF
                  IF(laieff_fit(ipts,ivm,ilevel)%b .LT. 0.0d0)THEN
                     ! temporary, during merge, to try and eliminate negative values
                     IF(printlev_loc.GT.4) &
                          WRITE(numout, *) '130416c laieff is at this point less than zero' 
                     laieff_fit(ipts,ivm,ilevel)%b = 0.0d0
                  END IF
                  IF(laieff_fit(ipts,ivm,ilevel)%c .LT. 0.0d0)THEN
                     ! temporary, during merge, to try and eliminate negative values
                     IF(printlev_loc.GT.4) &
                          WRITE(numout, *) '130416c laieff is at this point less than zero' 
                     laieff_fit(ipts,ivm,ilevel)%c = 0.0d0
                  END IF
                  IF(laieff_fit(ipts,ivm,ilevel)%d .LT. 0.0d0)THEN
                     ! temporary, during merge, to try and eliminate negative values
                     IF(printlev_loc.GT.4) &
                          WRITE(numout, *) '130416c laieff is at this point less than zero' 
                     laieff_fit(ipts,ivm,ilevel)%d = 0.0d0
                  END IF
                  IF(laieff_fit(ipts,ivm,ilevel)%e .LT. 0.0d0)THEN
                     ! temporary, during merge, to try and eliminate negative values
                     IF(printlev_loc.GT.4) &
                          WRITE(numout, *) '130416c laieff is at this point less than zero' 
                     laieff_fit(ipts,ivm,ilevel)%e = 0.0d0
                  END IF
          END DO
       ENDDO
    ENDDO

    ! Debug
    IF(printlev_loc.GT.4) THEN
       WRITE(numout, *) '130416b end of fitting_laieff'
       WRITE(numout, *) '130416b laieff_fit(1,4,:): ', laieff_fit(1,4,:)
    ENDIF
    !-

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving fitting_laieff'

  END SUBROUTINE fitting_laieff


!! ================================================================================================================================
!! SUBROUTINE 	: fit_laieff_to_points
!!
!>\BRIEF        
!!
!! DESCRIPTION  : 
!!
!!    NOTE:
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::laieff_fit
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE fit_laieff_to_points(nangles, test_angles, laieff, laieff_fit, ipts, ivm, ilevel)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(IN)                  :: nangles          !! The number of angles (data points)
    INTEGER(i_std),INTENT(IN)                  :: ipts             !! The grid square
    INTEGER(i_std),INTENT(IN)                  :: ivm              !! The PFT number
    INTEGER(i_std),INTENT(IN)                  :: ilevel           !! The canopy level
    REAL(r_std),DIMENSION(:),INTENT(IN)        :: test_angles      !! the cosine of the solar angles
    REAL(r_std),DIMENSION(:),INTENT(IN)        :: laieff           !! the laieff at each of these angles

    !! 0.2 Output variables
    TYPE(laieff_type),INTENT(out)              :: laieff_fit       !! Fitted parameters for the effective LAI

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)   :: x,y,one,x2,x3,x4,x5,x6,x7,x8,xy,x2y,x3y,x4y
    INTEGER  :: iangle
    REAL(r_std),DIMENSION(5) :: vector_b
    REAL(r_std),DIMENSION(5,5) :: matrix_a
    REAL(r_std) :: mean,sse,sst,r_corr, max_diff, diff

!_ ================================================================================================================================
    
    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering fit_laieff_to_points'


    ! We have a set of angles, test_angles.  We also have the
    ! effective LAI for each of these angles.  The structure
    ! laieff_fit gives the parameters and functional form.
    ! If you wish to change the functional form, you need to 
    ! change this routine, the lai_type, and the function
    ! below which evaluates the LAIeff as a function of
    ! the angle.
    
    ! We are going to use a least squares method to fit a parabola
    ! of the form a+b*x+c*x**2+d*x**3+e*x**4.  This means we will need to solve
    ! a system of linear equations, A x = b.
    ! In this notation, y is the sum over all angles for the laieff,
    ! x is the sum for the angles, x2 is the sum of the square of the
    ! angles, etc.

    y=zero
    xy=zero
    x2y=zero
    x3y=zero
    x4y=zero
    one=zero
    x=zero
    x2=zero
    x3=zero
    x4=zero
    x5=zero
    x6=zero
    x7=zero
    x8=zero
    DO iangle=1,nangles
           y=y+laieff(iangle)
           xy=xy+test_angles(iangle)*laieff(iangle)
           x2y=x2y+test_angles(iangle)**2*laieff(iangle)
           x3y=x3y+test_angles(iangle)**3*laieff(iangle)
           x4y=x4y+test_angles(iangle)**4*laieff(iangle)
           x=x+test_angles(iangle)
           x2=x2+test_angles(iangle)**2
           x3=x3+test_angles(iangle)**3
           x4=x4+test_angles(iangle)**4
           x5=x5+test_angles(iangle)**5
           x6=x6+test_angles(iangle)**6
           x7=x7+test_angles(iangle)**7
           x8=x8+test_angles(iangle)**8
    ENDDO
    one=REAL(nangles)

    ! Now we create our system of equations.
    vector_b(1:5) = (/ y, xy, x2y, x3y, x4y /)

    ! From what I see in the gauss_jordan routine, the format is
    ! matrix(row, column), which matches standard notation. In this
    ! case the matrix is symmetric, anyway.
    matrix_a(1,1:5) = (/ one,  x, x2, x3, x4 /) 
    matrix_a(2,1:5) = (/   x, x2, x3, x4, x5 /) 
    matrix_a(3,1:5) = (/  x2, x3, x4, x5, x6 /)
    matrix_a(4,1:5) = (/  x3, x4, x5, x6, x7 /)
    matrix_a(5,1:5) = (/  x4, x5, x6, x7, x8 /)

    ! TEST... should give -2, 3, 1...it does
!    matrix_a(1,1:3) = (/  -3, 2, -6 /) 
!    matrix_a(2,1:3) = (/  5, 7, -5 /) 
!    matrix_a(3,1:3) = (/  1, 4, -2 /)
!    vector_b(1:3) = (/ 6, 6, 8 /)

    CALL gauss_jordan_method(SIZE(vector_b),matrix_a,vector_b)

    ! These are our parameters.
    laieff_fit%a=vector_b(1)
    laieff_fit%b=vector_b(2)
    laieff_fit%c=vector_b(3)
    laieff_fit%d=vector_b(4)
    laieff_fit%e=vector_b(5)

    !++++++ CHECK ++++++++
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! These are some checks that we do to make sure our fitting
    ! is decent.  I'm going to leave these on during parameterization,
    ! but we might want to remove it during the real runs.  Things
    ! should be tested well before then, and it will only be freak
    ! chance that causes it to fail afterwards.  If this happens
    ! on one day in one hundred years, we don't care so much.

    ! Calculate a regression coefficienct
    sst=zero
    sse=zero
    mean=zero
    DO iangle=1,nangles
       mean=mean+laieff(iangle)
    ENDDO
    mean=mean/REAL(nangles)
    max_diff=zero
    DO iangle=1,nangles
       sst=sst+(laieff(iangle)-mean)**2
       diff=laieff(iangle)-calculate_laieff_fit(test_angles(iangle),laieff_fit)
       sse=sse+diff**2
       IF(diff .GT. max_diff) max_diff=diff
    ENDDO

    IF(sse .GT. sst)THEN
       r_corr = zero
    ELSEIF(sst .EQ. zero)THEN
       r_corr = zero
    ELSE
       r_corr=SQRT(un-sse/sst)
    ENDIF

    ! Debug
    IF(printlev_loc.GT.4 .AND. (ipts == test_grid) .AND. (ivm == test_pft)) &
         WRITE(numout,*) 'Regression coefficient: ',r_corr
    !-

    IF(r_corr .LT. 0.95)THEN

       ! Sometimes the correlation coefficient is not great, but the difference
       ! in the values is really small.  We don't care so much if our values
       ! are off by 0.001.

       ! Debug
       IF(printlev_loc.GT.4 .AND. (ipts == test_grid) .AND. (ivm == test_pft)) &
            WRITE(numout,*) 'max_diff ',max_diff
       !-

       IF(max_diff .GT. 0.01)THEN
          WRITE(numout,*) 'Difference is too big for fitting!'
          WRITE(numout,*) 'max_diff ',max_diff
          WRITE(numout,*) 'r_corr ',r_corr
          WRITE(numout,*) 'ipts,ivm,ilevel: ',ipts,ivm,ilevel
          WRITE(numout,*) 'COS(angle)     Real value       Fitted value'
          DO iangle=1,nangles
             WRITE(1600,'(10F20.10)') test_angles(iangle),&
                  laieff(iangle),calculate_laieff_fit(test_angles(iangle),laieff_fit)
             WRITE(numout,'(10F20.10)') test_angles(iangle),&
                  laieff(iangle),calculate_laieff_fit(test_angles(iangle),laieff_fit) 
          ENDDO
          CALL ipslerr_p (3,'fitting_laieff', 'Difference is too big for fitting.','','')
       ENDIF
    ENDIF
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Some debugging information
    IF(printlev_loc.GT.4 .AND. (ivm == test_pft) .AND. (ipts == test_grid))THEN
       ! Let's see how closely they match
       DO iangle=1,nangles
          IF(SUM(laieff(:)) .GT. 0.1 )THEN
             WRITE(numout,*) 'TESTING fitting LAIEFF ',ipts,ivm,ilevel
             WRITE(numout,'(10F16.8)') test_angles(iangle),&
                  laieff(iangle),calculate_laieff_fit(test_angles(iangle),laieff_fit)
          ENDIF
       ENDDO
    ENDIF

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving fit_laieff_to_points'

  END SUBROUTINE fit_laieff_to_points

!! ================================================================================================================================
!! SUBROUTINE 	: effective_lai
!!
!>\BRIEF        Calculate effective lai for the 1-D two stream 
!! radiative canopy transfer model
!!
!! DESCRIPTION  : The effective LAI is computed to be used with the two-stream
!!    albedo model described by Pinty 2006, by using an equation in Pinty 2004
!!    and the Pgap method of Haverd et al. 2012 code was taken from the latter group
!!    and cleaned up and reformatted...these routines were valided against the original
!!    code and gave identical results for Pgap (the original code was single precision, while
!!    these results are double precision):
!!
!!    NOTE:
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::laieff, ::Pgap_out
!!
!! REFERENCE(S) : Pinty et al; Simplifying the interaction of land surfaces 
!!   with radiation for relating remote sensing products to climate models
!!   JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D02116, doi:10.1029/2005JD005952, 2006
!!
!!                Pinty et al; Synergy between 1-D and 3-D radiation transfer models
!!   to retrieve vegetation canopy properties from remote sensing data
!!   JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 109, D21205, doi:10.1029/2004JD005214, 2004
!!
!!                Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!   (CanSPART) model: Formulation and application
!!   AGRICULTURAL AND FOREST METEOROLOGY, VOL. 160, pp. 14-35, 2012
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE effective_lai(npts, nlevels_loc, z_array_out,  circ_class_biomass, &
       circ_class_n, veget_max, coszang, lai_per_level, laieff, Pgap_cumul_out, Pgap_perlevel_out)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    INTEGER(i_std),INTENT(IN)                      :: nlevels_loc           !! The number of physical levels over which we
                                                                            !! will calculate the gap probability. We still 
                                                                            !! need  nlevels_loc in this subroutine call, 
                                                                            !! as this it also called from slowproc_veget, but
                                                                            !! only for the layer at the soil.
    INTEGER(i_std),INTENT(IN)                      :: npts                  !! Domain size - number of pixels (unitless)
                                                                            !! we will calculate the gap probability.
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(IN)      :: z_array_out           !! the height of the bottom of each of our levels @tex $(m)$ @endtex
                                                                            !! note that level 1 is closest to the ground
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN)  :: circ_class_biomass    !! Biomass of the componets of the model  
                                                                            !! tree within a circumference
                                                                            !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)      :: circ_class_n          !! Number of trees within each circ class
                                                                            !! class @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(IN)          :: veget_max             !! Maximum fraction of vegetation type including 
                                                                            !! non-biological fraction (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)          :: coszang               !! the cosine of the view zenith angle
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)      :: lai_per_level         !! LAI per vertical level
                                                                            !! @tex $(m^{2} m^{-2})$


    !! 0.2 Output variables

    REAL(r_std),DIMENSION(:,:,:),INTENT(OUT)       :: laieff                !! Leaf Area Index Effective converts 3D lai into 
                                                                            !! 1D lai for two stream radiation transfer model
                                                                            !! @tex $(m^{2} m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:,:),INTENT(OUT), OPTIONAL &
                                                   :: Pgap_cumul_out        !! The probability of finding a gap in the
                                                                            !! in canopy from the top of the canopy
                                                                            !! to a given level.
                                                                            !! (unitless)
    REAL(r_std),DIMENSION(:,:,:),INTENT(OUT), OPTIONAL &
                                                   :: Pgap_perlevel_out     !! The probability of finding a gap in the
                                                                            !! in the canopy between a level
                                                                            !! and the level immediately above.
                                                                            !! (unitless)


    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                   :: solar_angle           !! the view angle of every grid square (radians)
    INTEGER(i_std)                                 :: ilevel, ipts, icir    !! index (unitless)
    INTEGER(i_std)                                 :: ivm, jlevel           !! index (unitless)
    INTEGER(i_std),DIMENSION(npts,nvm)             :: upward_array          !! Flag to indicate if the view angle is pointing
                                                                            !! upwards (1) or downwards (0). (unitless)
 
    ! These next values need a bit of caution with the units.  Since they are derived
    ! from circ_class_biomass, it seems they have units per land area, in
    ! addition to the units given here.
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: crown_shadow_h        !! The projected shadow of an opaque crown.
                                                                            !! @tex $(m^2)$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: trunk_shadow_h        !! The projected shadow of the trunk.
                                                                            !! @tex $(m^2)$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: sbar_h                !! The mean free distance through the crown
                                                                            !! along the path dictated by the view angle.
                                                                            !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: favd_array            !! The foliage area volume density of
                                                                            !! the canopy.
                                                                            !! @tex $(LAI m^{-3})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: G_array               !! The mean horizontal projection of
                                                                            !! unit leaf area (leaf orientation function).
                                                                            !! @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,ncirc,nlevels_loc)      &
                                                   :: Pwc_h                 !! Crown porosity.
                                                                            !! (unitless)
    REAL(r_std),DIMENSION(npts,nvm,nlevels_loc)    :: TrunkShadowL          !! Trunk shadow, integrated over the population.
                                                                            !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,nlevels_loc)    :: Pbc_trunkL
    REAL(r_std),DIMENSION(nlevels_loc)             :: PorousCrownShadowL    !! The average shadow cast by a tree taking
                                                                            !! into account porosity.
                                                                            !! @tex $(m^2)$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,nlevels_loc)    :: PgapL_cumul           !! The probability of finding a gap in the
                                                                            !! in canopy from the top of the canopy
                                                                            !! to a given level.
                                                                            !! (unitless)
    REAL(r_std),DIMENSION(npts,nvm,nlevels_loc)    :: PgapL_perlevel        !! The probability of finding a gap in the
                                                                            !! in the canopy between a level
                                                                            !! and the level immediately above.
                                                                            !! (unitless)
    REAL(r_std),DIMENSION(npts,nvm)                :: favd                  !! The foliage area volume density of
                                                                            !! the canopy.
                                                                            !! @tex $(LAI m^{-3})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)               :: ind_loc               !! Density of individuals
                                                                            !! @tex $(m^{-2})$ @endtex
    ! The following two variables are ways to speed up the calculation
    ! by not trying to calculate everything, based on if we have
    ! biomass or individuals for this PFT on the grid square
    LOGICAL, DIMENSION(npts,nvm)                   :: lcalculate_tree       !! For trees
    LOGICAL, DIMENSION(npts,nvm)                   :: lcalculate_nontree    !! For non-trees

    REAL(r_std),DIMENSION(npts,nvm,ncirc,ndist_types)  :: values            !! An array that holds canopy properties.
    REAL(r_std)                                    :: lai_sum               !! The total LAI in the canopy,
                                                                            !! summed across all levels.
                                                                            !! @tex $(m^{2} m^{-2})$
    REAL(r_std),PARAMETER                          :: max_laieff=20.0       !! A numerical construct to cap the effective LAI so 
                                                                            !! as to not cause numerical problems. If laieff > 20
                                                                            !! the lai >> 20 so we have a problem somewhere else.
    REAL(r_std)                                    :: min_pgap              !! The probability that light makes
                                                                            !! it through a layer of max_laieff thickness.
    REAL(r_std),DIMENSION(npts,nvm)                :: favd_output           !! Variable to collect favd prior to writing it to the history files 
   
!_ ================================================================================================================================

    IF (firstcall_laieff) CALL ipslerr_p (3,'effective_lai', 'Initialization is not done before first call','','')

    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering effective_lai'

    !! 1. Initialize check for surface area conservation
    !  Veget_max is a INTENT(in) variable and can therefore
    !  not be changed during the course of this subroutine
    !  No need to check whether the subroutine preserves the
    !  total surface area of the pixel.

    !! 2. Initialize variables
    laieff(:,:,:) = zero
    PGapL_cumul(:,:,:) = un
    PGapL_perlevel(:,:,:) = un
    favd_output(:,:) = xios_default_val

    ! 0 means we are looking downward, which is what we need
    ! for the albedo since that is the direction the light is heading
    upward_array(:,:) = 0 

    ! Calculate the solar angle (radians)
    solar_angle(:)=ACOS(coszang(:)) 

    ! we need the total number of individuals
    ind_loc(:,:)=SUM(circ_class_n(:,:,:),3)
    ind_loc(:,1)=zero

    !! 3. Distinguish between tall and short vegetation
    ! The way Pgap is calculated differs between tall and short vegetation. Create
    ! a flag that can be used in subsequent calculations. For each PFT and each grid 
    ! square, let's determine if we actually need to run the calculation
    lcalculate_tree(:,:)=.FALSE.
    lcalculate_nontree(:,:)=.FALSE.
    DO ipts=1,npts
       DO ivm=1,nvm
          IF(is_tree(ivm))THEN
             IF(ind_loc(ipts,ivm) .GT. min_stomate .AND. &
                  ind_loc(ipts,ivm) .LT. val_exp*ncirc)THEN 
                lcalculate_tree(ipts,ivm)=.TRUE.
             ENDIF
          ELSE
             ! This approach makes the model crash if veget_max > 0 but there is
             ! no biomass yet. Replaced by an approach based on leaf biomass
             ! IF(veget_max(ipts,ivm) .GT. min_stomate) lcalculate_nontree(ipts,ivm)=.TRUE.
             IF(SUM(lai_per_level(ipts,ivm,:)) .GT. min_stomate) &
                 lcalculate_nontree(ipts,ivm)=.TRUE.

          ENDIF
       ENDDO
    ENDDO

    ! We need to derive some diagnostic variables from the biomass of each
    ! tree in the circ distribution, based on allometic relations. We need the
    ! height of the vegetation, the crown projected surface area, and the trunk
    ! diameter
    CALL derive_biomass_quantities(npts, nvm, ncirc, circ_class_n, &
         circ_class_biomass, values)

    !! 4. Calculate geometric quantities that go into the Pgap calculation
    DO ipts=1,npts
       DO ivm=1,nvm

          ! Calculate favd here such that we also calculate a value for the grass
          ! and croplands. For grass and croplands favd is NOT used in the calculation
          ! of Pgap but it is a useful variable to check the assumptions on the 
          ! crown shape and canopy height. Observed favd values for trees are below
          ! 3 m2/m3 but more often below 2m2/m3. Compute the average foliage volume 
          ! density across the height distribution. Notice that we do not need the 
          ! area of the grid square here because we have a number density of trees 
          ! instead of just the number of trees, so grid area cancels out in the equation
          CALL calculate_favd(ipts, ivm, values(ipts,ivm,:,icnvol), circ_class_biomass, &
               circ_class_n, favd)

          ! generate a larger array of values (we use the same favd for all layers and 
          ! diameter classes)
          favd_array(ipts,ivm,:,:) = favd(ipts,ivm)
          favd_output(ipts,ivm) = favd(ipts,ivm)

          IF (.NOT. lcalculate_tree(ipts,ivm)) CYCLE
          DO icir=1,ncirc
             
             ! These limits are arbitrary.
             IF(values(ipts,ivm,icir,icndiaver) .GT. laieff_zero_cutoff .AND. &
                  values(ipts,ivm,icir,icndiahor) .GT. laieff_zero_cutoff)THEN

                DO ilevel=1,nlevels_loc
               
                   ! compute the shadow of the crown projected onto each height for
                   ! each view angle using Appendix A in the manuscript...careful,
                   ! there are typos in the manuscript, but these routines came
                   ! directly from Vanessa and Jenny and should be correct if
                   ! either of the elipsoid axes are zero, this routine will crash.
                   ! Note that crown_shadow is defined for the canopy above the height
                   ! of the level (z_array_out). There is no crown_shadow above the
                   ! canopy.
                   crown_shadow_h(ipts,ivm,icir,ilevel) = spheroid_shadow(solar_angle(ipts),&
                        upward_array(ipts,ivm),&
                        values(ipts,ivm,icir,icndiaver),&
                        values(ipts,ivm,icir,icndiahor),&
                        values(ipts,ivm,icir,iheight),&
                        z_array_out(ipts,ivm,icir,ilevel))

                   ! Compute Eq. 10 in the manuscript, the mean path through
                   ! all canopy levels above z_array_out
                   sbar_h(ipts,ivm,icir,ilevel) = sbar(values(ipts,ivm,icir,iheight),&
                        values(ipts,ivm,icir,icndiaver)/2.0_r_std,&
                        values(ipts,ivm,icir,icndiahor)/2.0_r_std,&
                        z_array_out(ipts,ivm,icir,ilevel),&
                        crown_shadow_h(ipts,ivm,icir,ilevel),&
                        solar_angle(ipts), &
                        upward_array(ipts,ivm))

                END DO ! nlevels_loc

             ELSE

                ! Not enough LAI to bother with the calculations
                crown_shadow_h(ipts,ivm,icir,:)=zero
                sbar_h(ipts,ivm,icir,:) = zero
                
             ENDIF ! checking for values

!!$             ! This subroutine takes a lot of time, and we're not currently using 
!!$             ! the trunks below.  So let's just comment it out for now.
!!$             IF(values(ipts,ivm,icir,idiameter) .GT. laieff_zero_cutoff)THEN
!!$                
!!$                DO ilevel=1,nlevels_loc
!!$                
!!$                   ! Now compute the shadow of the trunks.  
!!$                   trunk_shadow_h(ipts,ivm,icir,ilevel) = trunk_shadow(solar_angle(ipts),&
!!$                        upward_array(ipts,ivm),&
!!$                        values(ipts,ivm,icir,idiameter),&
!!$                        values(ipts,ivm,icir,iheight),&
!!$                        z_array_out(ipts,ivm,icir,ilevel))
!!$
!!$                END DO ! loop over levels
!!$
!!$             ELSE
!!$
!!$                ! The trunk diameter is too small to throw a shadow
!!$                trunk_shadow_h(ipts,ivm,icir,ilevel)=zero
!!$
!!$             ENDIF

          ENDDO ! loop over circ classes 

!!$          ! do some integrations over the height distribution to find average values
!!$          ! This is equivalent to taking a weighted average by the number of 
!!$          ! individuals in each circ class.
!!$          DO ilevel=1,nlevels_loc
!!$             
!!$             ! find the average trunk shadow per square meter
!!$             ! this is Eq. 5 in the manuscript
!!$             TrunkShadowL(ipts,ivm,ilevel)=zero
!!$             DO icir=1,ncirc
!!$                ! Right now we are not using trunks.
!!$                TrunkShadowL(ipts,ivm,ilevel)=TrunkShadowL(ipts,ivm,ilevel)+&
!!$                     trunk_shadow_h(ipts,ivm,icir,ilevel)*&
!!$                     circ_class_n(ipts,ivm,icir)
!!$             ENDDO ! loop over circ classes
!!$          
!!$             ! this is part of Eq. 3 in the manuscript
!!$             TrunkShadowL(ipts,ivm,ilevel)=TrunkShadowL(ipts,ivm,ilevel)*&
!!$                  ind_loc(ipts,ivm)
!!$             Pbc_trunkL(ipts,ivm,ilevel)=EXP(-TrunkShadowL(ipts,ivm,ilevel))
!!$
!!$          ENDDO! loop over levels
    
          ! now we calculate the crown porosity, since the crown is not opaque...Eq. 8
          ! for this we need the projection of the leaf area in the direction of the beam,
          ! which is done by the Nilson method (end of Section 4.2 in the manuscript)
          ! the code also did this with the Ross G function, but we'll take Nilson's.
          ! the leaf angle distributions (LAD) are currently assumed to be spherical, which
          ! means the mean leaf inclination is 57° and the standard deviation is 20°,
          ! described with a two parameter beta function...this should be explored in more
          ! depth
          DO icir=1,ncirc

             DO ilevel=1,nlevels_loc
                ! this is done slightly differently from the original code, 
                ! because in the two-stream stream albedo method some coefficients
                ! are done under the assumption of Ross's G function and a
                ! spherical leaf distribution, which Pinty et al says makes the 
                ! G function constant and equal to 0.5. Using NilsonG and the
                ! above parameters gives a value of 0.505, so it's not much
                ! different but this is good for consistency. The original
                ! code reads: G_array(ipts,ivm,icir,ilevel)=&
                !  NilsonG(thetaL_mean_array(ipts,ivm,icir,ilevel),&
                !  thetaL_sd_array(ipts,ivm,icir,ilevel),laieff_solar_angle)
                G_array(ipts,ivm,icir,ilevel)=0.5_r_std

                ! I am adding a clumping factor here to account for needles 
                ! clumping to shoots, which increases the porosity of the canopy
                ! and lets more light through. Notice this is not in the 
                ! original Pgap model.
                Pwc_h(ipts,ivm,icir,ilevel)=EXP(-G_array(ipts,ivm,icir,ilevel)*&
                     favd_array(ipts,ivm,icir,ilevel)*sbar_h(ipts,ivm,icir,ilevel)&
                     *leaf_to_shoot_clumping(ivm))
                
                ! Error checking
                IF (Pwc_h(ipts,ivm,icir,ilevel).GT.un) THEN
                   WRITE(numout,*) 'ipts, ivm, icir, ilevel, ',ipts,ivm,icir,ilevel 
                   WRITE(numout,*) 'Pwc_h, ',Pwc_h(ipts,ivm,icir,ilevel)
                   WRITE(numout,*) 'G_array, ',G_array(ipts,ivm,icir,ilevel)
                   WRITE(numout,*) 'favd, ',favd_array(ipts,ivm,icir,ilevel)
                   WRITE(numout,*) 'sbar, ',sbar_h(ipts,ivm,icir,ilevel)
                   WRITE(numout,*) 'leaf_to_shoot_clumping, ',leaf_to_shoot_clumping(ivm)

                   ! hack_pwc_h can be used to accept this error and continue running
                   ! If hack_pwc_h is used, a small massbalance error will be seen in the next day
                   ! Therefor you will need to set ERR_ACT=1 to bypass this error
                   IF (hack_pwc_h) THEN
                      ! Send a warning and set PWc_h=1
                      WRITE(numout,*) 'Pwc_h>1 for ipts, ivm = ',ipts,ivm
                      CALL ipslerr_p(2,'stomate_laieff','Pwc_h .GT. 1','this is unexpected. Now set Pwc_h=1', &
                           "Note that this will lead to a small massbalance error in the following day, you'll need to set ERR_ACT=1")
                      Pwc_h(ipts,ivm,icir,ilevel)=un
                   ELSE
                      ! Stop the model
                      CALL ipslerr_p(3,'stomate_laieff','Pwc_h .GT. 1','this is unexpected','')
                   END IF
                END IF

             ENDDO ! loop over levels

          ENDDO ! loop over circ classes

       ENDDO ! loop over PFTs

    ENDDO ! loop over grid points

    ! Calculate the porous crown shadow integrating over the height distribution.
    ! This is equivalent to taking a weighted average by the number of individuals 
    ! in each circ class. Eq. 4, and then calculate Pgap
    DO ipts=1,npts
       DO ivm=1,nvm

          ! Here we need to calculate the effective LAI for all land types, 
          ! including grass, crops, and bare soil
          IF (is_tree(ivm)) THEN

             IF(.NOT. lcalculate_tree(ipts,ivm)) CYCLE

             ! first for forests
             DO ilevel=1,nlevels_loc
                PorousCrownShadowL(ilevel) = zero
                DO icir=1,ncirc
                   PorousCrownShadowL(ilevel) = PorousCrownShadowL(ilevel)+&
                        crown_shadow_h(ipts,ivm,icir,ilevel)*&
                        (1.0_r_std-Pwc_h(ipts,ivm,icir,ilevel))*&
                        circ_class_n(ipts,ivm,icir)
                ENDDO
                
                ! These next lines are commented out, but I'm leaving them here
                ! for clarity.  The first line normalises the integration done
                ! in the previous step so that we have the average porous
                ! crown shadow over the whole stand.  The second line multiplies
                ! by the density to be used in the exponential of the Poisson
                ! distribution.  It's clear from these two lines that
                ! they cancel each other out.
                ! PorousCrownShadowL(ipts,ivm,ilevel)=&
                !     PorousCrownShadowL(ipts,ivm,ilevel)/ind_loc(ipts,ivm)
                ! PorousCrownShadowL(ipts,ivm,ilevel)=&
                !     PorousCrownShadowL(ipts,ivm,ilevel)*ind_loc(ipts,ivm)
                
                ! This is the cumulative PgapL value. It represents the gap fraction
                ! considering all canopy above the level of interest (z_array_out)
                PGapL_cumul(ipts,ivm,ilevel) = EXP(-PorousCrownShadowL(ilevel))

                ! This line doesn't include the effect of trunks. Since our 
                ! trunks in the albedo calculation aren't distinguishable from
                ! the leaves, and since we don't want trunks absorbing light, we
                ! remove the influence of trunks for now.
                ! this seems to be Eq. 3 in the manuscript
!!$                PGapL_cumul(ipts,ivm,ilevel)=Pbc_trunkL(ipts,ivm,ilevel)*&
!!$                     EXP(-PorousCrownShadowL(ipts,ivm,ilevel))


             ENDDO ! loop over levels

          ELSEIF(ivm == ibare_sechiba)THEN

             ! For bare soil the transmission probability is always one. The
             ! only thing that matters to the albedo is the background reflectance
             PGapL_cumul(ipts,ivm,:)=un

          ELSE

             IF(.NOT. lcalculate_nontree(ipts,ivm)) CYCLE

             ! Now for grass and croplands. We treat these exactly the same 
             ! right now, just a slab of vegetation filling space.  For this, we
             ! need to know the amount of vegetation that is in each level, which
             ! was calculated in the find_lai_per_level routine.  It assumes folliage
             ! is distributed evenly from the ground to the height of the
             ! grass/crops.  Since the LAI will not be correct for grasslands and
             ! croplands due to the lack of management in ORCHIDEE, we are
             ! introducing a tunable parameter to help recover the true albedo values.

             DO ilevel=1,nlevels_loc
                ! This Eq. 1 in the manuscript, assuming G(theta)=0.5 (spherical
                ! distribution) and using the LAI value with the adjusted LAI 
                ! mentioned above. 
                ! Notice that this lai_temp value has to be cumulative, since that
                ! of the forests is and we correct Pgap based on that assumption below.
                ! This means that we need the total LAI through which light has passed
                ! in order to get here, i.e. the sum of all levels higher than and
                ! including this level.
                lai_sum=zero
                DO jlevel=nlevels_loc,ilevel,-1
                   lai_sum=lai_sum+lai_per_level(ipts,ivm,jlevel)
                ENDDO

                PGapL_cumul(ipts,ivm,ilevel)=EXP(-0.5_r_std*lai_sum*&
                     lai_correction_factor(ivm)/COS(solar_angle(ipts)))

                IF(ipts == test_grid .AND. ivm == test_pft .AND. printlev_loc.GT.4)THEN
                   WRITE(numout,*) 'PgapL cumulative: ',PGapL_cumul(ipts,ivm,ilevel)
                   WRITE(numout,*) 'lai_sum,lai_correction_factor(ivm): ',lai_sum,lai_correction_factor(ivm)
                   WRITE(numout,*) 'COS(solar_angle(ipts)),ilevel: ',COS(solar_angle(ipts)),ilevel
                ENDIF
             ENDDO

          ENDIF ! check for tree
       ENDDO ! loop over PFTs
    ENDDO ! loop over grid points

   
    ! NOTE: Find the pgap for each level seperately. This makes PGap a measure 
    ! of the chance to make it through just this level, as opposed to the 
    ! chance to make it through this level and all the above levels. 
    ! Converting Pgap from its initial cumulative value into a value representing
    ! a single layer rather than all layers above the height of interest is
    ! needed to calculate the effective lai for each layer. This effective
    ! lai is in turn used to calculate the albedo. Note that where we
    ! use Pgap to calculate veget, we need the cumalative Pgap at the
    ! bottom of the canopy.

    ! I do not need to CYCLE if lcalculate is not equal to
    ! true here because I initialize PgapL to be
    ! equal to 1 on all levels above.  The calculation
    ! of this loop is much less expensive than
    ! the other loops.
    DO ipts=1,npts
       DO ivm=1,nvm

          ! Debug
          IF(printlev_loc.GT.4)THEN
             IF((ipts == test_grid) .AND. (ivm == test_pft))THEN
                WRITE(numout,*) 'Cumulative Pgap'
                DO ilevel=nlevels_loc,1,-1
                   WRITE(numout,*) ivm,ilevel,PgapL_cumul(ipts,ivm,ilevel)
                ENDDO
             ENDIF
          ENDIF
          !-

          DO ilevel=1,nlevels_loc-1

             ! include a check for numerical stability.
             ! The minumum Pgap is a function of the maximum
             ! effective LAI that we allow.  The higher the
             ! effective LAI, the lower the amount of light
             ! which will be allowed through the layer.  This
             ! number cannot be so high as to cause numerical
             ! problems, though.  The value of 20 here is
             ! arbitrary, but doesn't seem to cause any problems.
             min_pgap=exp(-max_laieff/2.0/cos(solar_angle(ipts)))

             IF(PgapL_cumul(ipts,ivm,ilevel+1) .GT. min_pgap )THEN
                ! Calculate light that passes through an individual layer
                PgapL_perlevel(ipts,ivm,ilevel)= PgapL_cumul(ipts,ivm,ilevel)/ &
                     PgapL_cumul(ipts,ivm,ilevel+1)
             ELSE
                ! This is a problem.  We can't really say anything about
                ! the gap probability for this level.  However, if the
                ! above level is so dark, no light will make it
                ! to this level, anyway.  So let's just make this
                ! level opaque.  This is equivalent to giving a maximum
                ! effective LAI value.
                PgapL_perlevel(ipts,ivm,ilevel)=exp(-max_laieff/2.0/cos(solar_angle(ipts)))
             ENDIF

             ! Include a check.  PgapL should always be between zero and
             ! one.  Due to numerical issues, it's possible that it won't 
             ! be sometimes.  We don't check to see if it's less than
             ! zero because that will be caught later, and it's a more
             ! serious issue.
             IF(PgapL_perlevel(ipts,ivm,ilevel) .GT. un) PgapL_perlevel(ipts,ivm,ilevel)=un

          ENDDO ! loop over levels

          ! compute the effective LAI, based on Eq. 25 in Pintys 2004 paper cited
          ! in the references above
          ! Debug
          IF(test_pft == ivm .AND. test_grid == ipts .AND. printlev_loc.GT.4)THEN
             WRITE(numout,*) 'Effective lai',ivm
             WRITE(numout,*) 'ipts,ivm,nlevels_loc',ipts,ivm,nlevels_loc
          ENDIF
          !-

          DO ilevel=nlevels_loc,1,-1
             ! include a check for numerical stability...min_stomate is arbitrary
             IF(PgapL_perlevel(ipts,ivm,ilevel) .GT. min_stomate)THEN
                laieff(ilevel,ipts,ivm)=-2.0_r_std*COS(solar_angle(ipts))*&
                     LOG(PgapL_perlevel(ipts,ivm,ilevel))
             ELSE
                laieff(ilevel,ipts,ivm)=-2.0_r_std*COS(solar_angle(ipts))*&
                     LOG(min_stomate)
             ENDIF

             ! If the value is small enough or negative, we set it equal to zero.
             ! A negative value for the effective LAI does not make sense in the
             ! context of Eq. 12 in Pinty 2006.  If the effective LAI is zero, the
             ! uncollided transmission will be one.  If LAI is positive, the
             ! uncollided transmission will be between zero and one, since light
             ! is intercepted by the canopy elements.  A negative effective
             ! LAI results in uncollided transmitted light being greater than one,
             ! which means there is a light source in the canopy.
             IF(laieff(ilevel,ipts,ivm) .LE. min_stomate)THEN
                IF( laieff(ilevel,ipts,ivm) .LE. -min_stomate) THEN
                   WRITE(numout,*) 'ilevel,ipts,ivm', ilevel,ipts,ivm
                   WRITE(numout,*) 'laieff(ilevel,ipts,ivm)', laieff(ilevel,ipts,ivm)
                   CALL ipslerr_p (3,'effective_lai', 'We should not have a negative effective LAI.','','')
                ENDIF
                laieff(ilevel,ipts,ivm) = zero
             ENDIF


          ENDDO ! loop over levels
          
          ! Debug
          IF(test_pft == ivm .AND. test_grid == ipts .AND. printlev_loc.GT.4)THEN
             WRITE(numout,*) 'ilevel, laieff, PgapL (cumul), PgapL (per level): ',ivm
             DO ilevel=nlevels_loc,1,-1
                WRITE(numout,'(I5,3E20.10)') ilevel,laieff(ilevel,ipts,ivm),PgapL_perlevel(ipts,ivm,ilevel),PgapL_cumul(ipts,ivm,ilevel)
             ENDDO
          ENDIF
          !-

       ENDDO ! loop over PFTs
    ENDDO ! loop over grid points

    ! Write foliar area density to the history files
    CALL xios_orchidee_send_field("CROWN_POROSITY",favd_output(:,:))

    ! Deal with the optinal arguments
    IF(PRESENT(Pgap_cumul_out)) Pgap_cumul_out = PgapL_cumul
    IF(PRESENT(Pgap_perlevel_out)) Pgap_perlevel_out = PgapL_perlevel
    
    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving effective_lai'

  END SUBROUTINE effective_lai



!! ================================================================================================================================
!! FUNCTION 	: spheroid_shadow
!!
!>\BRIEF        Calculates the projected shadow of a tree canopy at different heights and
!!              sun angles, assuming a spheriod canopy
!!
!! DESCRIPTION  : Taken with permission from V. HAVERD's code, based on the equations in Appendix A
!!                of her paper
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::spheroid_shadow
!!
!! REFERENCE(S) : 
!!                Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  REAL(r_std) FUNCTION spheroid_shadow(theta,upward,T,D,&
       H,z)
    
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN)                      :: theta ! view angle (radians)
    REAL(r_std), INTENT(IN)                      :: T ! vertical diameter of the spheroid (m) 
    REAL(r_std), INTENT(IN)                      :: D ! horizontal diameter (m)
    REAL(r_std), INTENT(IN)                      :: H ! tree height (top of crown)
    REAL(r_std), INTENT(IN)                      :: z ! height above ground
    INTEGER(i_std), INTENT(IN)                   :: upward ! ==1 for upward view, ==0 for downward view
    
    
    !! 0.2 Output variables
!    REAL(r_std)                                  :: spheroid_shadow
    

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                  :: thetapr, Ypr, Y1pr, Y1, X1, D1pr, D2pr, Y

!_ ================================================================================================================================
    
    thetapr =  ATAN(T/D*TAN(theta)) ! zenith angle in spherical space
    
    
    ! the conventions are described in Appendix A of this paper
    IF (upward==1) THEN
       Y = z-H+T/2 
    ELSE
       Y = H-T/2-z
    ENDIF
    Ypr  = Y*D/T
    Y1pr = D/2.0*TAN(thetapr)/SQRT(1+TAN(thetapr)**2) 
    Y1 = T/D*Y1pr 
    IF (TAN(thetapr).NE.0.0) THEN
       X1 = Y1pr/TAN(thetapr) 
    ELSE
       X1 = D/2 
    END IF
    
    
    
    ! first calculate major axis (D1pr) of elliptic shadow
    IF (Y<=-T/2.0) THEN ! below crown bottom
        D1pr = 0.0 
    END IF
    
    IF (Y>-T/2 .AND. Y<=-Y1)  THEN !between crown bottom and lower tangent point
       D1pr = D*SQRT(1.0-Ypr**2*4.0/D**2)  
    END IF
    
    IF (Y>-Y1 .AND. Y <Y1) THEN !between lower and upper tangent points
       D1pr = D/2.0*SQRT(1.0-MIN(Ypr**2*4.0/D**2,1.0)) +X1 - &
            TAN(thetapr)*(-Y1pr-Ypr) 
    END IF
    
    IF (Y > Y1) THEN !above upper tangent point
       D1pr = D/COS(thetapr)     
    END IF

    ! second calculate minor axis (D2pr) of elliptic shadow 
    ! (this is just the diameter of the circular cross-section at height z)
    IF (Y <= -T/2) THEN
       D2pr = 0.0                      !shadow breadth
    END IF
    
    IF (Y<=0 .AND. Y >-T/2) THEN
       D2pr = 2.0*SQRT((D/2.0)**2-Ypr**2)    !shadow breadth
    END IF
    
    IF (Y>=0) THEN
       D2pr = D    !shadow breadth
    END IF
    

    ! shadow area
    spheroid_shadow = pi/4.0*D2pr*D1pr 

  END FUNCTION spheroid_shadow


!! ================================================================================================================================
!! FUNCTION     : trunk_shadow
!!
!>\BRIEF        ! Calculates the projected shadow of a tree trunk at different heights and
!!              sun angles, assuming a cone shape
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code, based on the equations in Appendix E
!!                of her paper
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : trunk_shadow
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION trunk_shadow(theta,upward,D,H,z)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN)                       :: theta ! view angle              @tex $(radians)$ @endtex
    REAL(r_std), INTENT(IN)                       :: D ! trunk diameter              @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)                       :: H ! tree height (top of crown)  @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)                       :: z ! height above ground         @tex $(m)$ @endtex
    INTEGER(i_std), INTENT(IN)                    :: upward ! ==1 for upward view, ==0 for downward view
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                  :: sinalpha, alpha 

!_ ================================================================================================================================
    
    sinalpha = D/(2.*H*TAN(theta))
    IF (ABS(sinalpha).LT.1.0.AND.D.GT.0.0) THEN
       alpha = 2.*ASIN(sinalpha)

       IF (upward==1) THEN
          ! shadow is bounded by base-of-trunk circle; trunk cross-section at z
          ! and tangents joining the two circles
          IF (z<H) THEN
             trunk_shadow = 1./4. * D**2*(1.-(1.-z/H)**2)/(TAN(alpha/2.)) + &
                  1./8.*D**2*(pi+alpha) + 1./8.*D**2*(1-z/H)**2*(pi-alpha)
          ELSE
             trunk_shadow = 1./4. * D**2*1./(TAN(alpha/2.)) + 1./8.*D**2*(pi+alpha)
          ENDIF
          
          trunk_shadow = trunk_shadow - pi * (D**2)/4.
       ELSEIF (upward==0) THEN
          
          ! shadow is bounded by trunk cross-section at z and tangents joining two circles
          IF (z<H) THEN
             trunk_shadow = 1./4. * D**2*(1.-z/H)**2*1./(TAN(alpha/2.)) + &
                  1./8.*D**2*(1.-z/H)**2*(pi+alpha)
          ELSE
             trunk_shadow = 0.
          ENDIF
          
       ENDIF
       
    ELSE
       trunk_shadow = 0.0
    ENDIF
    
  END FUNCTION trunk_shadow


!! ================================================================================================================================
!! FUNCTION     : sbar
!!
!>\BRIEF        ! Calculates the mean free distance through the crown, Eq. 10 in Haverd et al.
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : sbar
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION sbar(x, b, r, z, area,&
       theta, upward )
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN)         :: x             ! tree height              
                                                     ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)         :: b             ! vertical crown radius (diameter/2)
                                                     ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)         :: r             ! horizontal crown radius (diameter/2)
                                                     ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)         :: z             ! height above ground
                                                     ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)         :: theta         ! view angle
                                                     ! @tex $(radians)$ @endtex
    REAL(r_std), INTENT(IN)         :: area          ! shadow area of partial spheroid  
                                                     ! @tex $(m^2)$ @endtex
    INTEGER(i_std), INTENT(IN)      :: upward        ! ==1 for upward view, ==0 
                                                     ! for downward view
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                  :: Hpr, Y, V_ps

!_ ================================================================================================================================

    Hpr = x-b ! Hpr is height at crown centre
    IF (upward==1) THEN
       Y = z-Hpr  ! Y=0 corresponds to crown centre
    ELSE
       Y = Hpr - z
    ENDIF
    
    ! initialize the output variables, just in case there is a numerical reason 
    ! none of the loops is hit
    V_ps =zero
    sbar =zero

    IF (Y<b.AND.Y>-b) THEN

       ! volume of partial spheroid below height z [eq 69a in Haverd et al 2012]
       ! In haverd et al 2012 the equation for the volume of a 
       ! partial spheroid below height z is given. This volume has to
       ! be subtracted from the total volumd to obtain the volume
       ! above the height z.
       V_ps = 4.0_r_std/3.0_r_std*pi*(r)**2*(b) - ((2.0_r_std* b)/3.0_r_std - Y &
            + Y**3/(3.0_r_std* b**2)) *pi* r**2 

       ! there are some numerical issues here occasionally that give a negative 
       ! extremely small value, which in turn gives a negative sbar which causes
       ! problems later
       IF(ABS(V_ps) .LT. 0.0000001_r_std) V_ps= zero
       sbar = V_ps/COS(theta)/area

    ENDIF

   ! had to add the check for area [Eq. 69 in Haverd et al 2012]
   IF (Y>=b .AND. area .NE. zero) THEN
       V_ps = 4.0_r_std/3.0_r_std*pi*(r)**2*(b)  ; ! volume of whole spheroid
       sbar= V_ps/COS(theta)/area
    ENDIF

  END FUNCTION sbar
  
!! ================================================================================================================================
!! FUNCTION     : spheroid_area
!!
!>\BRIEF        ! Calculates the area of the cross section of a spheroid at a height of z
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : spheroid_area
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION spheroid_area(theta,upward,T,D,H,z)

    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN)          :: theta  ! view angle              
                                               ! @tex $(radians)$ @endtex
    REAL(r_std), INTENT(IN)          :: T      ! vertical diameter           
                                               ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)          :: D      ! horizontal diameter         
                                               ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)          :: H      ! tree height (top of crown)  
                                               ! @tex $(m)$ @endtex
    REAL(r_std), INTENT(IN)          :: z      ! height above ground         
                                               ! @tex $(m)$ @endtex
    INTEGER(i_std), INTENT(IN)       :: upward ! ==1 for upward view, 
                                               ! ==0 for downward view

  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                      :: thetapr, Ypr, Y1pr, Y1, X1, D2pr, Y

!_ ================================================================================================================================

    thetapr =  ATAN(T/D*TAN(theta)) ! zenith angle in spherical space

    IF (upward==1) THEN
       Y = z-H+T/2 
    ELSE
       Y = H-T/2-z
    ENDIF
    Ypr  = Y*D/T
    Y1pr = D/2.0*TAN(thetapr)/SQRT(1+TAN(thetapr)**2) 
    Y1 = T/D*Y1pr 
    IF (TAN(thetapr).NE.0.0) THEN
       X1 = Y1pr/TAN(thetapr) 
    ELSE
       X1 = D/2 
    END IF
    

    ! calculate minor axis (D2pr) of elliptic shadow 
    ! (this is just the diameter of the circular cross-section at height z)
    IF (Y<=T/2 .AND. Y >-T/2) THEN
       D2pr = 2.0*SQRT((D/2.0)**2-Ypr**2)    !shadow breadth
    ELSE
       D2pr = 0    !shadow breadth
    END IF

    ! shadow area
    spheroid_area = pi/4.0*D2pr**2 

  END FUNCTION spheroid_area

!! ================================================================================================================================
!! FUNCTION     : NilsonG
!!
!>\BRIEF        ! Calculates the projection of unit foliage area  G(theta) as described by
!!              Wang et al, 2007, referenced in section 4.2 from V. Haverd's paper
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : NilsonG
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION NilsonG(mean_thetaL,sd_thetaL,theta_deg)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN)          :: mean_thetaL ! the mean leaf angle inclination 
                                                    ! @tex $(degrees)$ @endtex 
    REAL(r_std), INTENT(IN)          :: sd_thetaL ! the standard deviation of the leaf angle inclination 
                                                  ! @tex $(degrees)$ @endtex 
    REAL(r_std), INTENT(IN)          :: theta_deg ! the solar zenith angle 
                                                  ! @tex $(degrees)$ @endtex 
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                      :: var_t, thl_bar, theta, thetaL, psi
    REAL(r_std), ALLOCATABLE         :: A(:), f(:)
    INTEGER(i_std)                   :: k, nthetaL
!_ ================================================================================================================================

    var_t = (sd_thetaL/90.0_r_std)**2
    thl_bar = mean_thetaL*pi/180.0_r_std
    theta = theta_deg*pi/180.0_r_std
    nthetaL = 20
    ALLOCATE (A(nthetaL))
    ALLOCATE (f(nthetaL))
    DO k=1,nthetaL
       thetaL = 0.0_r_std + 89.0_r_std/nthetaL*(k)*pi/180.0_r_std
       IF (ABS((1/TAN(theta))*(1/TAN(thetaL)))>1) THEN
          A(k) = COS(theta)*COS(thetaL);
       ELSE
          psi = ACOS((1.0_r_std/TAN(theta))*(1.0_r_std/TAN(thetaL)));
          A(k) = COS(theta)*COS(thetaL)*(1.0_r_std+2.0_r_std/pi*(TAN(psi)-psi));
       ENDIF
       f(k) = BETA_LAD(thl_bar,var_t, thetaL)
    ENDDO
    NilsonG = SUM(A(1:nthetaL)*f(1:nthetaL))/SUM(f(1:nthetaL))

    DEALLOCATE(A)
    DEALLOCATE(F)

  END FUNCTION NilsonG

!! ================================================================================================================================
!! FUNCTION     : BETA_LAD
!!
!>\BRIEF        ! Calculates the beta leaf angle distribution as described by
!!              Wang et al, 2007, referenced in section 4.2 from V. Haverd's paper
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : BETA_LAD
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION BETA_LAD(thl_bar,var_t, theta)

    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN):: thl_bar,var_t, theta
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                                  :: tbar, sig0_2, sigt_2, mu, nu, t, num
!_ ================================================================================================================================

    tbar = 2.0_r_std*thl_bar/(pi)
    sig0_2 = tbar*(1.0_r_std-tbar)
    sigt_2 = var_t
    num = sig0_2/sigt_2 - 1.0_r_std
    ! Test for negative num as this causes distribution to fail.
    IF (num <= zero) THEN 
       num = 0.01_r_std
       !write(*,*)'WARNING: leaf angle SD too large'
    ENDIF
    
    nu = tbar*num
    mu = (1.0_r_std-tbar)*num
    
    t = 2.0_r_std*theta/pi
    beta_lad = ((1.0_r_std-t)**(mu-1.0_r_std))*(t**(nu-1.0_r_std))/beta(mu,nu)   
    
  END FUNCTION BETA_LAD
!********************************************************************************

!! ================================================================================================================================
!! FUNCTION     : beta
!!
!>\BRIEF        ! Calculates the beta leaf angle distribution as described by
!!              Wang et al, 2007, referenced in section 4.2 from V. Haverd's paper
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : BETA
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION beta(z,w)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN):: z,w
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
!_ ================================================================================================================================

    beta=EXP(gammln(z)+gammln(w)-gammln(z+w))

  END FUNCTION beta
!! ================================================================================================================================
!! FUNCTION     : gammaln
!!
!>\BRIEF        ! Calculates the beta leaf angle distribution as described by
!!              Wang et al, 2007, referenced in section 4.2 from V. Haverd's paper
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : gammaln
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION gammln(xx)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std), INTENT(IN):: xx
  
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std) :: tmp,x,stp,temp,first,increment
    INTEGER(i_std) :: k,k2,n
    INTEGER(i_std), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    REAL(r_std),DIMENSION(6) :: coef,ARTH
!_ ================================================================================================================================
    
    stp = 2.5066282746310005_r_std
    coef = (/76.18009172947146_r_std,&
         -86.50532032941677_r_std,24.01409824083091_r_std,&
         -1.231739572450155_r_std,0.1208650973866179e-2_r_std,&
         -0.5395239384953e-5_r_std/)
    x=xx
    tmp=x+5.5_r_std
    tmp=(x+0.5_r_std)*LOG(tmp)-tmp
    
    n = SIZE(coef)
    increment = 1.0_r_std
    first = x+1.0_r_std
    IF (n > 0) arth(1)=first
    IF (n <= NPAR_ARTH) THEN
       DO k=2,n
          ARTH(k)=ARTH(k-1)+increment
       END DO
    ELSE
       DO k=2,NPAR2_ARTH
          ARTH(k)=ARTH(k-1)+increment
       END DO
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       DO
          IF (k >= n) EXIT
          k2=k+k
          arth(k+1:MIN(k2,n))=temp+arth(1:MIN(k,n-k))
          temp=temp+temp
          k=k2
       END DO
    END IF
    
    gammln=tmp+LOG(stp*(1.000000000190015_r_std+&
         SUM(coef(:)/arth(:)))/x)
  END FUNCTION gammln

!! ================================================================================================================================
!! FUNCTION     : partial_spheroid_vol
!!
!>\BRIEF        ! Calculates the beta leaf angle distribution as described by
!!              Wang et al, 2007, referenced in section 4.2 from V. Haverd's paper
!!
!! DESCRIPTION : Taken with permission from V. HAVERD's code
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : partial_spheroid_vol
!!
!! REFERENCE(S) : Haverd et al; The Canopy Semi-analytic Pgap and radiative transfer 
!!      (CanSPART) model: Formulation and application, AGRICULTURAL AND FOREST METEOROLOGY
!!      Vol. 160, pp. 14-35, 2012
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION partial_spheroid_vol(Z_up, Z_down, b, r, h)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    ! The units of all these variables can be anything, so long as they
    ! are all the same.
    !  I assume
    REAL(r_std), INTENT(IN)         :: Z_up  ! The height of the top cutting plane
    REAL(r_std), INTENT(IN)         :: Z_down ! The height of the lower cutting plane
    REAL(r_std), INTENT(IN)         :: b     ! Distance from center to top of the spheroid
    REAL(r_std), INTENT(IN)         :: r     ! Distance from center to outside on the XY plane
    REAL(r_std), INTENT(IN)         :: h     ! Height of the top of the spheroid above the ground
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std)                      :: Z1, Z2, center

!_ ================================================================================================================================
    
    ! initialize the output variables
    partial_spheroid_vol =zero

    ! It's possible that none of this spheroid is in this level!
    center = h-b

    IF(Z_down - center .GE. b)RETURN
    IF(Z_up - center .LE. -b)RETURN

    ! Z1 and Z2 have to be within the limits of the spheroid
    ! I also want the center at the origin.
    Z2=MIN(Z_up-center,b)
    Z1=MAX(Z_down-center,-b)

    ! With a little calculus, I computed the following formula, which
    ! gives the correct result in the case of the full spheroid,
    ! i.e. Z2=b and Z1=-b, as well as for a full sphere (r=b)
    partial_spheroid_vol = pi*r**2*(Z2-Z1-(Z2**3-Z1**3)/(3.0_r_std*b**2))

    ! there are some numerical issues here occasionally that give a negative 
    ! extremely small value, which in turn gives a negative partial_spheroid_vol which causes
    ! problems later
    IF(ABS(partial_spheroid_vol) .LT. min_stomate) partial_spheroid_vol= zero

  END FUNCTION partial_spheroid_vol

!! ================================================================================================================================
!! SUBROUTINE 	: calculate_favd
!!
!>\BRIEF       Calculates the foliage area volume density of a canopy.
!!
!! DESCRIPTION  : WARNING: This should only be used for trees.  Grasses
!!                and crops have no canopy volume currently in the model.
!!
!!    NOTE:
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::favd
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE calculate_favd(ipts, ivm, cn_vol, circ_class_biomass, &
       circ_class_n, favd)

  !! 0 Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER(i_std),INTENT(in)                            :: ipts         
    INTEGER(i_std),INTENT(in)                            :: ivm        
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)        :: circ_class_biomass    !! Biomass of the components of the model  
                                                                                  !! tree within a circumference class
                                                                                  !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)            :: circ_class_n          !! Number of trees within each circ class
                                                                                  !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)                :: cn_vol                !! Canopy volume for each circ class
                                                                                  !! @tex $(m^{-3} ind^{-1})$ @endtex
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: favd                  !! The LAI per canopy volume density
                                                                                  !! @tex $(m^{-3} m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER :: icir
    REAL(r_std) :: total_bm
    REAL(r_std) :: total_cn_vol
    REAL(r_std) :: max_height
    REAL(r_std) :: total_la
    REAL(r_std) :: lai_temp
    REAL(r_std) :: ind_loc
!_ ================================================================================================================================
    
    IF (printlev_loc >= 4) WRITE(numout,*) 'Entering calculate_favd'

    ind_loc = SUM(circ_class_n(ipts,ivm,:))

    IF(is_tree(ivm))THEN
       ! How much total leaf biomass do we have?
       ! How much total canopy volume do we have?
       total_bm=zero
       total_cn_vol=zero
       max_height=zero
       
       ! Maximum height
       max_height = MAXVAL(wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
            ivm,pipe_tune2(ipts,ivm)))

       ! Notice that we don't have to take the surface area of the stand
       ! into account.  total_bm will be in units of gC / m**2, while
       ! the canopy volume should be in units of m**3 / m**2.  The surface
       ! area would cancel out when we take the ratio.
       DO icir=1,ncirc

          ! Right now, total_cn_vol is in m**3/tree * tree/m**2
          total_cn_vol=total_cn_vol+cn_vol(icir)*&
               circ_class_n(ipts,ivm,icir)
          !            circ_class_n(ipts,ivm,icir)*surface_area

          ! Debug
          IF(ipts == test_grid .AND. ivm == test_pft .AND. printlev_loc.GT.4)THEN
             WRITE(numout,'(A,I5,3F16.8)') 'icir,cn_vol,circ_class_biomass,circ_class_n: ',&
                  icir,cn_vol(icir),circ_class_biomass(ipts,ivm,icir,ileaf,icarbon),&
                  circ_class_n(ipts,ivm,icir)
          ENDIF
          !-

       ENDDO

       ! total_cn_vol is calculated in m3/m2 so its unit is m (similar to how 
       ! precipitation is expressed). If total_cn_vol > max_height we clearly 
       ! have a problem.
       total_cn_vol = MIN(total_cn_vol,max_height)

       ! How much leaf area do we have on our grid square? Make use
       ! of cc_to_lai to correctly account for sla dynamics.
       ! Total_la is unitless. m**2 / m**2.
       total_la = cc_to_lai(circ_class_biomass(ipts,ivm,:,ileaf,icarbon),&
            circ_class_n(ipts,ivm,:),ivm)
       
       ! Debug
       IF(ipts == test_grid .AND. ivm == test_pft .AND. printlev_loc.GT.4)THEN
          WRITE(numout,'(A,3F16.8)') &
               'total_la,total_cn_vol: ',total_la,total_cn_vol 
       ENDIF
       !-

       ! What is the foliage area volume density?
       IF (total_cn_vol .GT. min_stomate) THEN
          ! 1 / m**3 (canopy volume) / m**2 (land area)
         favd(ipts,ivm)=total_la/total_cn_vol
       ELSE
          favd(ipts,ivm)=zero
       ENDIF

    ELSE

       ! What is the height of our grass/crop?
       lai_temp = biomass_to_lai(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)*circ_class_n(ipts,ivm,:)), ivm)

       ! This converts LAI to a height.  I took this from functional_allocation.
       max_height = lai_temp * lai_to_height(ivm)

       ! The units here are per square meter, so the volume of the cube is
       ! just the height*1*1.
       IF (max_height .GT. min_stomate) THEN
          favd(ipts,ivm)=lai_temp / max_height
       ELSE
          favd(ipts,ivm)=zero
       ENDIF

    ENDIF

    IF (printlev_loc >= 5) WRITE(numout,*) 'Leaving calculate_favd'

  END SUBROUTINE calculate_favd

!! ================================================================================================================================
!! FUNCTION     : calculate_laieff_fit
!!
!>\BRIEF        : 
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : calculate_laieff_fit
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  REAL(r_std) FUNCTION calculate_laieff_fit(test_angle, laieff_fit)
    !! 0 Variable and parameter declaration 
  
    !! 0.1 Input variables
    REAL(r_std),INTENT(IN)        :: test_angle        !! the cosine of the solar angle
    TYPE(laieff_type),INTENT(in)    :: laieff_fit      !! Fitted parameters for the effective LAI

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
!_ ================================================================================================================================

    IF (firstcall_laieff) CALL ipslerr_p (3,'calculate_laieff_fit', 'Initialization is not done before first call','','')

    ! y = a + b*x + c*x**2 + d*x**3 + e*x**4
    calculate_laieff_fit=laieff_fit%a + test_angle * laieff_fit%b + &
         test_angle**2 * laieff_fit%c +  test_angle**3 * laieff_fit%d + &
         test_angle**4 * laieff_fit%e



  END FUNCTION calculate_laieff_fit

END MODULE stomate_laieff
