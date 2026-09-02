! =================================================================================================================================
! MODULE       : sapiens_forestry
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Gathers the main elements for forest management: the "sapiens_forestry_main"
!! subroutine, which itself calls a set of subroutines
!! sapiens_forestry_clear, clearcut, thinning, harvest, force_load, QsortC, and
!! Partition, and a set of functions used in these subroutines.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)   :
!! - Asael, S., 1999. Typologie des peuplements forestiers du massif vosgiens. 
!! C.R.P.F. Lorraine-Alsace, Nancy, 54 p.
!! - Bellassen, V., Le Maire, G., Dhote, J.F., Viovy, N., Ciais, P., 2010. 
!! Modeling forest management within a global vegetation model – Part 1: 
!! model structure aSnd general behaviour. Ecological Modelling 221, 2458–2474.
!! - Bellassen, V., Le Maire, G., Guin, O., Dhote, J.F., Viovy, N., Ciais, P., 
!! 2011a. Modeling forest management within a global vegetation model – Part 2: 
!! model validation from tree to continental scale. Ecological Modelling 222, 
!! 57–75.
!! - Bellassen, V., Viovy, N., Luyssaert, S., Le Maire, G., Schelhaas, M.J., 
!! Ciais, P., 2011b. Reconstruction and attribution of the carbon sink of 
!! European forests between 1950 and 2000. Global Change Biology 17, 3274-3292.
!! - Bottcher, H., Kurz, W.A., Freibauer, A., 2008. Accounting of forest carbon
!! sinks and sources under a future climate protocol-factoring out past 
!! disturbance and management effects on age-class structure. Environmental 
!! Science & Policy 11, 669-686.
!! - Cazin, A., Vallet, P., Dhote, J.F., 2003. Propositions LERFoB pour CARBOFOR
!! Atelier modélisation.
!! - Condes, S., Sterba, H., 2005. Derivation of compatible crown width 
!! equations for some important tree species of Spain. Forest Ecology and 
!! Management 217, 203-218.
!! - Deleuze, C., Pain, O., Dhote, J.F., Herve, J.C., 2004. A flexible radial
!! increment model for individual trees in pure even-aged stands. Annals of 
!! Forest Science 61, 327-335.
!! - Dhôte, J.-F., Hervé, J.-C., 2000. Changements de productivité dans quatre 
!! forêts de chênes sessiles depuis 1930 : une approche au niveau du peuplement. 
!! Ann. For. Sci. 57, 651-680.
!! - Dhôte, J.F., 1999. Compétition entre classes sociales chez le chêne sessile
!! et le hêtre. Revue Forestière Française, 309-325.
!! - Dhôte, J.F., Le Moguédec, G., 2003. Présentation du modèle Fagacées. 
!! INRA, 31 p.
!! - IFN, 2008. Raw inventory data - 2005-2007 inventory campaigns, www.ifn.fr.
!! - Kolari, P., Pumpanen, J., Rannik, U., Ilvesniemi, H., Hari, P., Berninger,
!! F., 2004. Carbon balance of different aged Scots pine forests in Southern 
!! Finland. Global Change Biology 10, 1106-1119.
!! - Lanier, L., 1994. Précis de sylviculture. Ecole Nationale du Génie Rural, 
!! des Eaux et des Forêts (ENGREF), Nancy, 477 p.
!! - Law, B.E., Sun, O.J., Campbell, J., Van Tuyl, S., Thornton, P.E., 2003. 
!! Changes in carbon storage and fluxes in a chronosequence of ponderosa pine. 
!! Global Change Biology 9, 510-524.
!! - Liberloo, M., Calfapietra, C., Lukac, M., Godbold, D., Luos, Z.B., Polle, 
!! A., Hoosbeek, M.R., Kull, O., Marek, M., Raines, C., Rubino, M., Taylor, G., 
!! Scarascia-Mugnozza, G., Ceulemans, R., 2006. Woody biomass production during 
!! the second rotation of a bio-energy Populus plantation increases in a future 
!! high CO2 world. Global Change Biology 12, 1094-1106.
!! - Liberloo, M., Luyssaert, S., Bellassen, V., Djomo, S.N., Lukac, M., 
!! Calfapietra, C., Janssens, I., Hoosbeek, M.R., Viovy, N., Churkina, G., 
!! Scarascia-Mugnozza, G., Ceulemans, R., 2010. Bio-energy retains its 
!! mitigation potential under elevated CO2. PLOS-one 5, e1 1648.
!! - Litton, C.M., Raich, J.W., Ryan, M.G., 2007. Carbon allocation in forest 
!! ecosystems. Global Change Biology 13, 2089-2109.
!! - Luyssaert, S., Inglima, I., Jung, M., Richardson, A.D., Reichsteins, M., 
!! Papale, D., Piao, S.L., Schulzes, E.D., Wingate, L., Matteucci, G., Aragao, 
!! L., Aubinet, M., Beers, C., Bernhoffer, C., Black, K.G., Bonal, D., 
!! Bonnefond, J.M., Chambers, J., Ciais, P., Cook, B., Davis, K.J., Dolman, 
!! A.J., Gielen, B., Goulden, M., Grace, J., Granier, A., Grelle, A., Griffis,
!! T., Grunwald, T., Guidolotti, G., Hanson, P.J., Harding, R., Hollinger, D.
!! Y., Hutyra, L.R., Kolar, P., Kruijt, B., Kutsch, W., Lagergren, F., 
!! Laurila, T., Law, B.E., Le Maire, G., Lindroth, A., Loustau, D., Malhi, Y.,
!! Mateus, J., Migliavacca, M., Misson, L., Montagnani, L., Moncrieff, J., 
!! Moors, E., Munger, J.W., Nikinmaa, E., Ollinger, S.V., Pita, G., 
!! Rebmann, C., Roupsard, O., Saigusa, N., Sanz, M.J., Seufert, G., Sierra, C., 
!! Smith, M.L., Tang, J., Valentini, R., Vesala, T., Janssens, I.A., 2007. CO2 
!! balance of boreal, temperate, and tropical forests derived from a global 
!! database. Global Change Biology 13, 2509-2537.
!! - Mokany, K., Raison, R.J., Prokushkin, A.S., 2006. Critical analysis of 
!! root: shoot ratios in terrestrial biomes. Global Change Biology 12, 84-96.
!! - Mund, M., Kummetz, E., Hein, M., Bauer, G.A., Schulze, E.D., 2002. Growth
!! and carbon stocks of a spruce forest chronosequence in central Europe. 
!! Forest Ecology and Management 171, 275-296.
!! - Newton, R.F., Amponsah, I.G., 2007. Comparative evaluation of five height-
!! diameter models developed for black spruce and jack pine stand-types in 
!! terms of goodness-of-fit, lack-of-fit and predictive ability. Forest Ecology
!! and Management 247, 149-166.
!! - Numerical Recipes, 2007. Numerical Recipes: The Art of Scientific Computing. 
!! Cambridge University Press, 1256 p.
!! - Ovington, J.D., Madgwick, H.A.I., 1957. Afforestation and soil reaction. 
!! Journal of Soil Science 8, 141-149.
!! - Pontailler, J.Y., Ceulemans, R., Guittet, J., 1999. Biomass yield of poplar 
!! after five 2-year coppice rotations. Forestry 72, 157-163.
!! - Pretzsch, H., Biber, P., Durský, J., 2002. The single tree-based stand 
!! simulator SILVA: construction, application and evaluation. Forest Ecology and 
!! Management 162, 3.
!! - Reineke, L.H., 1933. Perfecting a stand-density index for even-aged forests. 
!! Journal of Agricultural Research 46, 627-638.
!! - Turner, D.P., Ritts, W.D., Cohen, W.B., Maeirsperger, T.K., Gower, S.T., 
!! Kirschbaum, A.A., Running, S.W., Zhao, M.S., Wofsy, S.C., Dunn, A.L., Law, B.E.,
!! Campbell, J.L., Oechel, W.C., Kwon, H.J., Meyers, T.P., Small, E.E., Kurc, S.A.,
!! Gamon, J.A., 2005. Site-level evaluation of satellite-based global terrestrial 
!! gross primary production and net primary production monitoring. Global Change 
!! Biology 11, 666-684.
!! - Vacchiano, G., Motta, R., Long, J.N., Shaw, J.D., 2008. A density management 
!! diagram for Scots pine (Pinus sylvestris L.): A tool for assessing the forest's 
!! protective effect. Forest Ecology and Management 255, 2542-2554.
!! - Vieira, I.C.G., de Almeida, A.S., Davidson, E.A., Stone, T.A., de Carvalho, 
!! C.J.R., Guerrero, J.B., 2003. Classifying successional forests using Landsat 
!! spectral properties and ecological characteristics in eastern Amazonia. Remote 
!! Sensing of Environment 87, 470-481.
!! - Zianis, D., Mencuccini, M., 2004. On simplifying allometric analyses of forest 
!! biomass. Forest Ecology and Management 187, 311-332.
!! 
!! SVN      :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/sapiens_forestry.f90 $
!! $Date: 2026-01-29 14:53:46 +0100 (jeu. 29 janv. 2026) $
!! $Revision: 9352 $
!! \n
!_ ================================================================================================================================

MODULE sapiens_forestry

  ! modules used:
#ifdef CPP_IEEE_ARITHMETIC
  USE,INTRINSIC :: IEEE_ARITHMETIC, only : IEEE_IS_NAN
#endif

  USE netcdf
  USE ioipsl_para
  USE constantes
  USE grid
  USE pft_parameters
  USE function_library,     ONLY: wood_to_stand_volume, wood_to_height, &
                                  wood_to_circ, Nmax, &
                                  wood_to_dia, sort_circ_class_biomass, &
                                  wood_to_qmdia, wood_to_qmdia_up_half, calculate_rdi, &
                                  check_vegetation_area, check_mass_balance, get_printlev, &
                                  calculate_rdi_boundaries
  USE interpol_help,        ONLY: aggregate_p
  USE stomate_data
  USE xios_orchidee
  USE time,                 ONLY: FirstDayYear
  
  IMPLICIT NONE

  ! private & public routines

  PRIVATE

  PUBLIC sapiens_forestry_main,                sapiens_forestry_clear,               &
         sapiens_forestry_set_fm,              sapiens_forestry_set_species_change, &
         sapiens_forestry_set_desired_fm,      &
         sapiens_forestry_litter_raking,       sapiens_forestry_species_change,      &
         sapiens_forestry_flag_species_change, sapiens_forestry_xios_initialize
  
  LOGICAL, SAVE                 :: firstcall_sapiens_forestry = .TRUE.   !! first call
!$OMP THREADPRIVATE(firstcall_sapiens_forestry)  
  LOGICAL, SAVE                 :: xios_interpolation_fm                 !! Flag for interpolating using XIOS 
!$OMP THREADPRIVATE(xios_interpolation_fm)
  INTEGER(i_std), SAVE          :: printlev_loc                          !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
  CHARACTER(LEN=80), SAVE       :: filename_fm                           !! Filename read from run.def for fm file
!$OMP THREADPRIVATE(filename_fm)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_clear
!!
!>\BRIEF        Set the flag ::firstcall_sapiens_forestry to .TRUE. and as such activate section 2 of the subroutine forestry (see below).
!! 
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_clear
    firstcall_sapiens_forestry = .TRUE.
  END SUBROUTINE sapiens_forestry_clear


  !!  =============================================================================================================================
  !! SUBROUTINE:    sapiens_forestry_xios_initialize
  !!
  !>\BRIEF	  Initialize xios dependant defintion before closing context defintion
  !!
  !! DESCRIPTION:	  Initialize xios dependant defintion needed for the interpolations done in sapiens_forestry_read_fm
  !!                      Reading is deactivated if the stomate restart file exists because the variable
  !!                      should be in the restart file already.
  !!                      This subruting is called before closing context with xios_orchidee_close_definition in intersurf 
  !!                      via the subroutine sechiba_xios_initialize. 
  !!
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE sapiens_forestry_xios_initialize
    
    CHARACTER(LEN=255)   :: name            !! Filename without suffix .nc
    LOGICAL              :: lerr            !! Flag to dectect error
     
 
    !Config Key   = FM_FILE
    !Config Desc  = Name of file from which the forest management map is to be read
    !Config If    = OK_READ_FM_MAP
    !Config Def   = FMmap.nc
    !Config Help  = The name of the file to be opened to read a forest management
    !Config         map (including a layer for every PFT) is given here. 
    !Config Units = [FILE]    filename_fm = 'FMmap.nc'
    filename_fm = 'FMmap.nc'
    CALL getin_p('FM_FILE',filename_fm)

    ! Remove suffix .nc from filename and define new name to XIOS
    name = filename_fm(1:LEN_TRIM(filename_fm)-3)
    CALL xios_orchidee_set_file_attr("fm_file",name=name)


    !Config Key   = XIOS_INTERPOLATION_FM
    !Config Desc  = Read and interpolate using XIOS
    !Config If    = OK_READ_FM_MAP
    !Config Def   = XIOS_INTERPOLATION
    !Config Help  = This key is proposed for debugging purposes when only one file will be read by XIOS
    !Config Units = [FILE]
    xios_interpolation_fm = xios_interpolation
    CALL getin_p('XIOS_INTERPOLATION_FM',xios_interpolation_fm)
    
    ! Check if the forest management file will be read by XIOS, by IOIPSL or not at all
    IF ( (xios_interpolation_fm .AND. ok_read_fm_map) .AND. &
         ( (stom_restname_in=='NONE' .OR. vegetmap_reset) .OR. &
            (FirstDayYear .AND. veget_update>0) ) ) THEN
       ! The forest management file will be red using XIOS
       IF (printlev>=1) WRITE(numout,*) 'Reading of forest managemend file will be done later using XIOS. The filename is ', filename_fm
    ELSE
       ! The file will be red with the method using IOIPSL or not at all
       ! Some prints for clairification
       IF ( .NOT. ok_read_fm_map) THEN 
          IF (printlev>=1) WRITE (numout,*) 'The forest management file will not be read becuase OK_READ_FM_MAP=n'
       ELSE IF ((stom_restname_in=='NONE' .OR. vegetmap_reset) .OR. &
            (FirstDayYear .AND. veget_update>0) ) THEN
          IF (printlev>=1) WRITE (numout,*) 'The forest management file will be red later using IOIPSL (not using XIOS)'
       ELSE
          IF (printlev>=1) WRITE (numout,*) "OK_READ_FM_MAP=y, the forest management file will be red from restart file"
       END IF

       ! The file will not be read by XIOS. Now deactivate the reading for XIOS.
       CALL xios_orchidee_set_file_attr("fm_file",enabled=.FALSE.)
       CALL xios_orchidee_set_fieldgroup_attr("fm_frac",enabled=.FALSE.)
    END IF
    
  END SUBROUTINE sapiens_forestry_xios_initialize

!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_main
!!
!>\BRIEF        This subroutine is the core of the forest management module. It
!! has been modified so it only decides if trees are getting cut or not.
!!
!! DESCRIPTION   : This forest management module has been changed quite a bit
!! from Valentin's original version (documented in Ballassen et al (2010).  The
!! biggest reason for these changes is that it is meant to be used with the new
!! functional allocation scheme.  This scheme already includes circumference
!! classes.  In each circumference class there are circ_class_n model trees
!! which are all identical.  Therefore, we don't need to create a picture
!! of every tree on the stand since in reality we only have ncirc trees, each
!! of which are replicated.
!!
!! In addition to this, several other processes that were done here were moved
!! to other parts of the code.  For example, trees are just scheduled to die
!! here; the actually killing of the biomass and moving the wood to the harvest
!! and litter pools is done in sapiens_kill.f90.  All tree establishment after
!! clearcutting is done in stomate_prescribe.f90.  Self-thinning is a completely
!! natural process, and therefore done with environmental mortality in 
!! stomate_mark_kill.f90.

!!++++ CHECK REFERENCE ++++
!! Liberloo et al. (2010) describes how short 
!! rotation coppices are simulated and validated.\n
!!+++++++++++++++++++++++

!! The forest management flags have changed as well.
!! FM = 1 : No human intervention (ORCHIDEE default)
!!      2 : Thinnings based on the RDI, clearcuts based on tree density,
!!              annual increment, and tree diameter.  Thinnings from above and
!!              from below are determined by the sign of thstrat.
!!      3 : Coppices
!!      4 : Short rotation coppices
!!      5 : Continuous cover forestry
!!
!! Relative density index (RDI) target of human thinning at the beginning of 
!! the rotation (unitless). This is the value aimed at by human thinning at 
!! the beginning of the rotation. It is necessarily lower than 1. The actual 
!! RDI varies around this target, + or - delta_rdi, when forest_managed = 2. 
!! When forest_managed = 3, human thinning occurs every year.  For FM = 4,
!! the thinnings occur as a function of stand age.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_kill
!!
!! REFERENCE(S)   : See above, module description.
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_main (npts, age_stand, age_stand_bm, last_cut, &
             circ_class_n, circ_class_kill, forest_managed, &
             circ_class_biomass, mai, pai, previous_wood_volume, &
             mai_count, coppice_dens, veget_max, fm_change_map, &
             species_change_map, whychange, tot_trees)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: npts                !! Domain size - number of pixels 
                                                                                  !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max           !! "maximal" coverage fraction of a PFT on 
                                                                                  !! the ground (unitless, 0-1) 
    INTEGER(i_std), DIMENSION(:,:),INTENT(in)              :: last_cut            !! Years since last thinning (years)
    INTEGER(i_std), DIMENSION(:,:),INTENT(in)              :: age_stand           !! Age of stand (years)
    REAL(r_std), DIMENSION(:,:),INTENT(in)                 :: age_stand_bm        !! Conserved biomass-weighted stand age (yr) - STAND_AGE, DIAGNOSTIC only since 2026-07-27
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)             :: fm_change_map       !! A map which gives the desired FM strategy when
                                                                                  !! the PFT will be replanted after a clearcut.
                                                                                  !! (1-nvm,unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)             :: species_change_map  !! A map which gives the PFT number that each
                                                                                  !! PFT will be replanted as in case of a clearcut.
                                                                                  !! (1-nvm,unitless)
    !! 0.2 Output

    !! 0.3 Modified fields
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)       :: circ_class_biomass  !! Biomass components of the model tree  
                                                                                  !! within a circumference class
                                                                                  !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: circ_class_n        !! Number of trees in each circumference 
                                                                                  !! class (0-2.2+ and poles-large wood, see 
                                                                                  !! part 12.) 
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)          :: forest_managed      !! Forest management flag: 1= Unmanaged (ORCHIDEE default),
                                                                                  !! 2= Thinnings, 3= Coppices, 4= Short rotation
                                                                                  !! coppices, 5= Continuous cover forestry
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)       :: circ_class_kill     !! Number of trees within a circ that needs
                                                                                  !! to be killed @tex $(ind m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: mai                 !! The mean annual increment
                                                                                  !! @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: pai                 !! The period annual increment
                                                                                  !! @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: previous_wood_volume!! The volume of the tree trunks
                                                                                  !! in a stand for the previous year.
                                                                                  !! @tex $(m**3 / m**2 )$ @endtex
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)           :: mai_count           !! The number of times we've
                                                                                  !! calculated the volume increment
                                                                                  !! for a stand
    REAL(r_std), DIMENSION(:,:),INTENT(inout)              :: coppice_dens        !! The density of a coppice at the first
                                                                                  !! cutting.
                                                                                  !! @tex $( 1 / m**2 )$ @endtex

    REAL(r_std), DIMENSION(:,:),INTENT(inout)              :: whychange          
    REAL(r_std), DIMENSION(:,:),INTENT(inout)              :: tot_trees

    !! 0.4 Local variables

    INTEGER(i_std)                                         :: icir,ivm,ipts,ifm   !! Indexes
    INTEGER(i_std)                                         :: imbc, iele, ipar    !! Indexes
    REAL(r_std)                                            :: trees,total_trees   !! a number of trees
    REAL(r_std), DIMENSION(ncirc)                          :: diameters_temp      !! Trunk diameter calculated from allometry [m]
    REAL(r_std), DIMENSION(ncirc)                          :: circ_temp           !! Trunk circumference calculated from allometry
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                 :: height,dia,cn       !! Tree height calculated from allometric 
                                                                                  !! relationships (m)
    REAL(r_std)                                            :: ave_tree_height
    REAL(r_std)                                            :: woody_biomass
    REAL(r_std)                                            :: current_wood_volume
    REAL(r_std)                                            :: current_increment
    REAL(r_std)                                            :: ave_tree_dia
    INTEGER(i_std)                                         :: nrotations          !! The number of rotations that we've done
                                                                                  !! for SRC minus one
    REAL(r_std)                                            :: qm_dia              !! quadratic mean diameter of the forest (m)
    REAL(r_std)                                            :: qm_dia_up_half      !! quadratic mean diameter of the dominant class (m)
    REAL(r_std)                                            :: rot_target_loc      !! AGE_ROTATION: effective per-pixel target rotation age (yr): map value or scalar
    REAL(r_std), DIMENSION(npts)                           :: ph_f_cut            !! PROGRESSIVE_HARVEST: rotation-paced fraction of the mature area to cut (-)
    REAL(r_std), DIMENSION(npts,nvm)                       :: ph_f_cut_pft        !! HARVEST_FEED_COUPLING: per-PFT cut fraction (group quota); flag off = pixel value copied to every PFT (-)
    REAL(r_std)                                            :: ph_a_grp            !! HARVEST_FEED_COUPLING: thinning-managed area of the species group (m2)
    REAL(r_std)                                            :: ph_a_term           !! HARVEST_FEED_COUPLING: area of the group terminal age class (m2)
    REAL(r_std)                                            :: ph_r4               !! HARVEST_FEED_COUPLING: terminal-class filling ratio vs MAT_TARGET_FRAC (-)
    REAL(r_std)                                            :: ph_gmod             !! HARVEST_FEED_COUPLING: demand modulation g, in [harv_g_min, harv_g_max] (-)
    INTEGER(i_std)                                         :: ph_jv, ph_jv1, ph_jv2 !! HARVEST_FEED_COUPLING: group scan index and PFT bounds of the group
    REAL(r_std)                                            :: ph_a_forest         !! PROGRESSIVE_HARVEST: forest area of the pixel (m2)
    REAL(r_std)                                            :: ph_demand           !! PROGRESSIVE_HARVEST: surface a couper par an, Sigma_jv A_jv/R(jv) (m2)
    REAL(r_std)                                            :: ph_a_mature         !! PROGRESSIVE_HARVEST: area of the slots meeting the clearcut criterion (m2)
    REAL(r_std)                                            :: ph_slot_area        !! PROGRESSIVE_HARVEST: ground area of one PFT slot (m2)
    REAL(r_std)                                            :: ph_qm_dia_up_half   !! PROGRESSIVE_HARVEST: dominant-class QMD in the pre-pass (m)
    LOGICAL                                                :: ph_active           !! PROGRESSIVE_HARVEST: this slot is cut fractionally this year
    LOGICAL                                                :: ph_dia_ok           !! PROGRESSIVE_HARVEST: eligibilite de coupe (test de diametre retire sous PH)
    LOGICAL                                                :: thin_this_slot      !! FALSE dans la derniere classe d age : phase de recolte, pas d eclaircie
    LOGICAL                                                :: forced_thin_slot    !! TRUE si le creneau subit l eclaircie ANNUELLE forcee
    INTEGER(i_std)                                         :: iage_ft             !! Rang de la classe de maturite dans son groupe, pour l eclaircie forcee
    REAL(r_std)                                            :: rdi_thin_goal       !! RDI vise par l eclaircie : bande basse, ou cible forcee
    INTEGER(i_std)                                         :: iage_th             !! Indice de classe d age pour la bande de RDI fixe
    INTEGER(i_std)                                         :: iage_thin_lo        !! Premiere classe d age eclaircie : 2, ou 1 sous OK_RDI_BAND_CLASS1
    REAL(r_std)                                            :: dfac_loc            !! Modulation du diametre de coupe par l'intensite de gestion (-)
    LOGICAL                                                :: allow_mgmt_cut      !! AGE_CLASS_BOUNDS_PFT: this slot may undergo a MANAGEMENT clearcut (last age class)
    REAL(r_std), PARAMETER                                 :: ph_min_survivor = 0.01 !! PROGRESSIVE_HARVEST: minimum surviving share below which the cut is treated as a whole-slot clearcut (-)
    INTEGER(i_std)                                         :: ierr                !! error message
    INTEGER(i_std)                                         :: igroup_new          !! Species group to which the species with which
                                                                                  !! the stand will be replanted belongs (unitless)
    INTEGER(i_std)                                         :: igroup_old          !! Species group to which the current species
                                                                                  !! belongs (unitless)
    REAL(r_std), DIMENSION(ncirc)                          :: circ_class_n_new    !! variable used to sort circ_class_n
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)         :: circ_class_biomass_new !! variable used to sort circ_class_biomass
    REAL(r_std)                                            :: ramp_u              !! Avancement du peuplement dans sa classe d'amelioration (0-1)
    REAL(r_std)                                            :: ramp_lo             !! Consigne de RDI interpolee sur la rampe (-)
    LOGICAL                                                :: cond1_collapsed     !! CONDITION 1 : le peuplement est-il EFFONDRE (et non ouvert a dessein)
    REAL(r_std), DIMENSION(npts,nvm)                       :: rdi                 !! Relative density index (unitless, 0-2)
    REAL(r_std), DIMENSION(npts,nvm)                       :: rdi_target_upper    !! Upper limit of RDI. When reached the stand
                                                                                  !! needs to be thinned (unitless)   
    REAL(r_std), DIMENSION(npts,nvm)                       :: rdi_target_lower    !! Lower limit of RDI. When thinned, thin to this             
                                                                                  !! stand density (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)     :: check_intern        !! Contains the components of the internal
                                                                                  !! mass balance chech for this routine
                                                                                  !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)             :: closure_intern      !! Check closure of internal mass balance
                                                                                  !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)             :: pool_start          !! Start and end pool of this routine 
                                                                                  !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)             :: pool_end            !! Start and end pool of this routine 
                                                                                  !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    INTEGER(i_std), DIMENSION(npts,nvm)                    :: spinup_clearcut     !! Map to indicate clearcut event during spinup 
                                                                                  !! for a given PFT and pixel (zero = no clearcut; one = clearcut).
    INTEGER(i_std), DIMENSION(npts,nvmap)                  :: fm_map_temp         !! A temporary variable to hold the forest
    INTEGER(i_std)                                         :: spinup_clearcut_force
!_ ===============================================================================================================================

 !! 1. Initialize
 
   IF (firstcall_sapiens_forestry) THEN   
       ! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_forestry=.FALSE.
    END IF

    IF ( printlev_loc >= 2) WRITE(numout,*) 'Entering sapiens_forestry_main'

    !! 1.2 Initialize check for mass balance closure
    !  We don't actually have to do this in sapiens_forestry_main, since none of the 
    !  pools are changed even touched. The sole purpose of this routine 
    !  is to schedule trees for killing which means that only 
    !  ::circ_class_kill is changed. Biomass is sorted so we will check
    !  wether nothing goes wrong during sorting
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements

          DO ipar = 1,nparts

             DO icir = 1,ncirc

                ! Initial biomass pool. Use circ_class_biomass
                ! because that variable is the prognostic variable
                ! biomass is syncronized to circ_class_biomass. If 
                ! many diameter classes are used, each diameter class
                ! separatly passes all criteria but when combined
                ! a mass balance problem occurs.
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                   circ_class_biomass(:,:,icir,ipar,iele) * &
                   circ_class_n(:,:,icir) * veget_max(:,:)

             ENDDO
             
          ENDDO

       ENDDO

    ENDIF ! err_act.GT.1

    !! 1.3 Initialize check for surface area conservation
    !  Veget_max is a INTENT(in) variable and can therefore
    !  not be changed during the course of this subroutine
    !  No need to check whether the subroutine preserves the
    !  total surface area of the pixel.

    !! 1.4 Sort the diameter classes
    ! A small number of individuals and distributing small labile and 
    ! carbres pools from the stand level (at which they are calculated) to 
    ! the tree level (at which they are stored) can result in precision
    ! errors. Sometimes these errors are sufficient to overturn the ascending
    ! order of tree biomass that is expected in circ_class_biomass. 
    ! sapiens_forestry_main is called after this routine and relies on the 
    ! assumption that the biomass is increasing. When needed, sort biomass
    ! to make sure this assumption is satisfied. Sorting is enforced at
    ! several places in the code
    DO ipts = 1,npts
       DO ivm = 1,nvm
          IF (is_tree(ivm)) THEN
             circ_class_biomass_new = circ_class_biomass(ipts,ivm,:,:,:)
             circ_class_n_new = circ_class_n(ipts,ivm,:)
             CALL sort_circ_class_biomass(circ_class_biomass_new, &
                  circ_class_n_new)
             circ_class_biomass(ipts,ivm,:,:,:) = circ_class_biomass_new
             circ_class_n(ipts,ivm,:) = circ_class_n_new
          ENDIF
       ENDDO
    ENDDO
    

    !! 1.5 Initialize spinup_clearcut from file or run.def
    !!     Note: this is done in sapiens_forestry_main because it is only called once a year
    ! We need to have the option to read the forest management
    ! strategy from a map (NetCDF file).  If this option is
    ! equal to Y, we will overwrite the forest_managed_forced 
    ! option above, so you should be careful to only use one
    ! or the other.
    IF (ok_read_sp_clearcut_map) THEN
       ! If we are using age classes, we read in the map in the same
       ! way but then we change it a bit to account for age classes.
       
       CALL sapiens_forestry_read_spinup_clearcut(npts, lalo, neighbours, resolution, &
            contfrac, fm_map_temp)

       IF (nagec .GT. 1)THEN
          ! All age classes of the same PFT will go through clearcut
          ! at the same time during spinup
          DO ivm = 1,nvm
             spinup_clearcut(:,ivm)=fm_map_temp(:,agec_group(ivm))
          ENDDO
       ELSE
          spinup_clearcut(:,:)=fm_map_temp(:,:)
       ENDIF

    ELSE
       
       ! Read value from run.def
       
       !Config Key   = SPINUP_CLEARCUT_FORCE
       !Config Desc  = New species after a final cut for testing and debugging only
       !Config If    = .NOT. OK_READ_SP_CLEARCUT_MAP
       !Config Def   = 0
       !Config Help  = 
       !Config Units = [0/1]
       spinup_clearcut_force=0
       CALL getin_p('SPINUP_CLEARCUT_FORCE',spinup_clearcut_force)
       
       spinup_clearcut(:,:)=spinup_clearcut_force
    ENDIF
   
    ! Output spinup_clearcut
    CALL xios_orchidee_send_field("SPINUP_CLEARCUT",REAL(spinup_clearcut))

    !! 1.6 PROGRESSIVE_HARVEST: rotation-paced harvest cadence
    !  Guillaume M. -- A forest on rotation R renews 1/R of its area per year, a constraint
    !  the all-or-nothing clearcut below ignores: it empties whichever slot meets the
    !  trigger. The pace comes from R, not from the clearcut size S, which only sets the
    !  texture on the ground. See design/MODULE_DESIGN_PROGRESSIVE_HARVEST.md section 2.2.
    !
    !      f_cut(ji) = MIN[ 1 , ( A_forest(ji) / R(ji) ) / A_mature(ji) ]
    !
    !  Guillaume M. -- A_mature is the area of the slots satisfying the CONDITION-2 trigger,
    !  so it must be known BEFORE any cut is applied -> this pre-pass. MIN[1,.] degrades
    !  gracefully: when the mature pool is smaller than the annual due we take all that is
    !  available and carry no debt forward.
    !
    !  Guillaume M. -- /!\ Deliberate approximation: the pre-pass replicates CONDITION 2 ONLY.
    !  Slots exiting earlier through CONDITION 1/1bis stay in A_mature, which is a
    !  denominator setting a pace, not a cut list.
    IF (ok_progressive_harvest) THEN

       IF (.NOT. ALLOCATED(ph_cut_frac)) THEN
          ALLOCATE(ph_cut_frac(npts,nvm), stat=ierr)
          IF (ierr /= 0) CALL ipslerr_p(3,'sapiens_forestry_main', &
               'Pb allocation ph_cut_frac','','')
       ENDIF

       ! Guillaume M. -- Reset every year: ph_cut_frac describes THIS year's cut and is
       ! consumed downstream by sapiens_kill (harvest_area scaling) then age_class_distr
       ! (slot split).
       ph_cut_frac(:,:) = zero
       ph_f_cut(:) = zero
       ph_f_cut_pft(:,:) = zero

       IF (ok_harvest_feed_coupling) THEN

          ! Guillaume M. -- HARVEST_FEED_COUPLING (design doc 4.4): one quota PER species
          ! group (working circle), demand and supply both filtered on ifm_thin, demand
          ! modulated by the terminal-class filling ratio through the MAT_BAND dead band.
          ! The historical pixel-level quota below is kept verbatim as the flag-off path.
          DO ipts = 1,npts
             DO ph_jv = 1,nvm

                IF (.NOT. is_tree(ph_jv)) CYCLE
                ! One pass per group, entered at its first PFT
                IF (ph_jv .NE. start_index(agec_group(ph_jv))) CYCLE
                ph_jv1 = ph_jv
                ph_jv2 = ph_jv + nagec_pft(agec_group(ph_jv)) - 1

                ph_demand   = zero
                ph_a_mature = zero
                ph_a_grp    = zero
                ph_a_term   = zero
                DO ivm = ph_jv1, ph_jv2
                   IF (veget_max(ipts,ivm) == zero) CYCLE
                   ! Guillaume M. -- Demand and supply share the SAME management filter:
                   ! the historical path counts every tree PFT in the demand but only
                   ! ifm_thin mature slots in the supply, inflating f_cut wherever
                   ! management is mixed. Fixed here only -- flag off keeps it as is.
                   IF (forest_managed(ipts,ivm) .NE. ifm_thin) CYCLE
                   ph_slot_area = veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
                   ph_a_grp = ph_a_grp + ph_slot_area
                   IF (ALLOCATED(target_rotation_age)) THEN
                      ph_demand = ph_demand + ph_slot_area &
                           / MAX(target_rotation_age(ipts,ivm), min_stomate)
                   ENDIF
                   ! r4 counts the terminal-class AREA (an empty mature slot is still area)
                   IF (ivm .EQ. ph_jv2) ph_a_term = ph_a_term + ph_slot_area
                   ! Supply keeps the main-loop guards: no biomass, nothing to cut
                   IF (SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .LE. min_stomate) CYCLE
                   IF (ok_clearcut_last_class .AND. ivm .NE. ph_jv2) CYCLE
                   ph_a_mature = ph_a_mature + ph_slot_area
                ENDDO ! ivm

                IF (ph_a_mature .LE. min_stomate .OR. ph_demand .LE. min_stomate) CYCLE

                ! Guillaume M. -- Dead band shared with the maturity regulator: inside it
                ! the rotation rules alone. A single-class group has no maturity structure
                ! to key on -> g stays 1. gain = 0 gives g = 1 everywhere (granularity-only
                ! mode); the MAX(r4, min_stomate) keeps 0**gain well defined.
                ph_gmod = un
                IF (nagec_pft(agec_group(ph_jv)) .GT. 1 .AND. ph_a_grp .GT. min_stomate) THEN
                   ph_r4 = ph_a_term / (MAX(mat_target_frac(nagec_pft(agec_group(ph_jv))), &
                        min_stomate) * ph_a_grp)
                   IF (ph_r4 .LT. mat_band_low .OR. ph_r4 .GT. mat_band_high) THEN
                      ph_gmod = MIN(harv_g_max, MAX(harv_g_min, &
                           MAX(ph_r4, min_stomate)**harv_coupling_gain))
                   ENDIF
                ENDIF

                ph_f_cut_pft(ipts, ph_jv1:ph_jv2) = MIN(un, ph_gmod * ph_demand / ph_a_mature)

             ENDDO ! ph_jv
          ENDDO ! ipts

       ELSE

       DO ipts = 1,npts

          ph_a_forest = zero
          ph_a_mature = zero
          ph_demand   = zero

          DO ivm = 1,nvm

             IF (.NOT. is_tree(ivm)) CYCLE
             IF (veget_max(ipts,ivm) == zero) CYCLE

             ph_slot_area = veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
             ph_a_forest = ph_a_forest + ph_slot_area
             ! Guillaume M. -- Harvest demand is accumulated PER PFT: the rotation is per
             ! PFT, so the area to cut is Sum_jv A_jv/R(jv) and not A_forest/R.
             ! /!\ This must stay INSIDE the ivm loop: read after the ENDDO,
             ! target_rotation_age would be indexed with a stale ivm.
             IF (ALLOCATED(target_rotation_age)) THEN
                ph_demand = ph_demand + ph_slot_area &
                     / MAX(target_rotation_age(ipts,ivm), min_stomate)
             ENDIF

             ! Guillaume M. -- Same guards as the main loop: no biomass, nothing to manage.
             IF (SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .LE. min_stomate) CYCLE

             ! Guillaume M. -- CONDITION 2 applies to ifm_thin only (see the trigger below).
             IF (forest_managed(ipts,ivm) .NE. ifm_thin) CYCLE

             ! Guillaume M. -- The management-cut gate must apply here too, otherwise
             ! A_mature would count slots the main loop then refuses to cut and f_cut
             ! would be under-estimated by the same amount.
             IF (ok_clearcut_last_class) THEN
                IF (ivm .NE. start_index(agec_group(ivm)) &
                     + nagec_pft(agec_group(ivm)) - 1) CYCLE
             ENDIF

             ! Guillaume M. -- Being in the last age class is enough to allow the cut: that
             ! class carries the "mature stand liable to be harvested" function. The cut
             ! diameter still acts, but in a single place -- the lower bound of that class,
             ! age_class_bound(nagec-1): crossing the bound means entering the class, and
             ! entering the class means becoming cuttable.
             ph_a_mature = ph_a_mature + ph_slot_area

          ENDDO ! ivm

          IF (ph_a_mature .GT. min_stomate .AND. ALLOCATED(target_rotation_age)) THEN

             IF (ph_demand .GT. min_stomate) THEN

                ! Guillaume M. -- MIN(un, ...) is the ONLY bound: f_cut is an area fraction.
                ! PH_MAX_FRACTION / PH_MIN_FRACTION were removed as no-ops -- the cap
                ! duplicated this MIN and the floor tested f_cut < 0 on a ratio of positive
                ! areas. See constantes.f90.
                ph_f_cut(ipts) = MIN(un, ph_demand / ph_a_mature)

             ENDIF

          ENDIF

          ! The consumption sites read ph_f_cut_pft: copying the pixel value to every
          ! PFT keeps the flag-off path bit-identical to the historical scalar.
          ph_f_cut_pft(ipts,:) = ph_f_cut(ipts)

       ENDDO ! ipts

       ENDIF ! ok_harvest_feed_coupling

    ENDIF ! ok_progressive_harvest


 !! 2. Check and apply sapiens_forestry_main
 
    ! for every forest stand in the model, we need to decide 
    ! if we are thinning or clearcutting.
    DO ipts=1,npts

       pft:DO ivm=1,nvm

          ! Forest management can only be applied to forests
          IF(.NOT. is_tree(ivm)) CYCLE

          ! This PFT is not present, so it can't be managed
          IF(veget_max(ipts,ivm) == zero) CYCLE

          ! There is no biomass, so nothing to manage
          IF(SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .LE. min_stomate) CYCLE

          ! Guillaume M. -- Under ok_clearcut_last_class, MANAGEMENT clearcuts (CONDITION
          ! 2/2bis and 3/3bis) are restricted to the last age class, which carries the
          ! end-of-rotation function. Deliberately not covered: technical resets, the
          ! stand-failure CONDITION 1, species change on coppice and every natural
          ! disturbance -- all must be able to reset any class.
          allow_mgmt_cut = .TRUE.
          IF (ok_clearcut_last_class) THEN
             allow_mgmt_cut = ( ivm .EQ. start_index(agec_group(ivm)) &
                  + nagec_pft(agec_group(ivm)) - 1 )
          ENDIF
          
          ! For unmanaged forests, either we apply clearcut during spinup, or 
          ! simply do nothing.
          IF ( forest_managed(ipts,ivm) .EQ. ifm_none ) THEN

             IF ( (spinup_clearcut(ipts,ivm) .EQ. flag_spinup_clearcut) .OR. forced_clear_cut ) THEN
                CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                    circ_class_kill(ipts,ivm,:,ifm_none,icut_clear))
             ELSE
                ! Otherwise we just get out of the loop
                CYCLE
             ENDIF

          ELSE

             ! Managed pft, apply to appropriate management
             IF (printlev_loc>=4) THEN
                WRITE(numout,*) 'Sapiens_forestry_main: check whether the stand ',&
                     'needs management', ivm
                CALL flush(numout)
             ENDIF

             ! Convert the aboveground woody biomass into a wood volume, 
             ! removing the influence of the branches.
             current_wood_volume = wood_to_stand_volume(npts,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:),ivm,branch_ratio(ivm),0)

             ! The mean annual increment is tricky. Need to calculate it
             ! for the volume, since that is what foresters use and the NPP
             ! would include leaf NPP.  However, we cannot calculate a
             ! volume increment immediately after thinning.  Instead, the first
             ! year after a thinning or after a clearcut (for consistency)
             ! we just save the current volume.
             IF(last_cut(ipts,ivm) .EQ. 1)THEN

                ! do nothing.  We'll save the current wood volume later.
                
             ELSEIF(last_cut(ipts,ivm) .GT. 1)THEN

                current_increment=&
                     current_wood_volume-previous_wood_volume(ipts,ivm)
                mai_count(ipts,ivm)=mai_count(ipts,ivm)+1
                mai(ipts,ivm)=mai(ipts,ivm)*REAL(mai_count(ipts,ivm)-1)/&
                     REAL(mai_count(ipts,ivm))+current_increment/&
                     REAL(mai_count(ipts,ivm))

                ! The periodic annual increment is an average over the past few
                ! years.  Notice that the first year after a cut that we can
                ! compute an increment for is year 2 (the difference between 
                ! year 2 and year 1).     
                IF(last_cut(ipts,ivm) .LT. (n_pai+1))THEN
                   
                   pai(ipts,ivm)=pai(ipts,ivm)*REAL(last_cut(ipts,ivm)-2)/&
                        REAL(last_cut(ipts,ivm)-1)+&
                        current_increment/REAL(last_cut(ipts,ivm)-1)

                ELSE

                   pai(ipts,ivm)=pai(ipts,ivm)*REAL(n_pai-1)/REAL(n_pai)+&
                        (current_wood_volume-previous_wood_volume(ipts,ivm))&
                        /REAL(n_pai)

                ENDIF

             ELSE

                ! Cannot think of a scenario where last_cut is zero.
                ! If it is zero, we thinned or clearcut, but that happens
                ! after this routine and we always increment it before this
                ! routine.
                WRITE(numout,*) 'ERROR: Sapiens_forestry_main, why is last_cut equal to zero?'
                WRITE(numout,*) 'ipts, ivm, last_cut(ipts,ivm) ',&
                     ipts, ivm, last_cut(ipts,ivm)
                IF (err_act.GT.1) CALL ipslerr_p (3, &
                     'sapiens_forestry_main', 'last_cut is zero.',&
                     'This should not happen.',&
                     'Look in the output file for ERROR.')
             ENDIF

             ! We also need to store the wood volume of this year in order to 
             ! calculate the PAI for next year.
             previous_wood_volume(ipts,ivm)=current_wood_volume

             ! our rdi_target will become a function of both PFT 
             ! and diameter/density Calculate rdi and rdi_target 
             ! and delta_rdi based on the yield models rdi of the 
             ! stand as it is
             qm_dia = wood_to_qmdia(&
                     circ_class_biomass(ipts,ivm,:,:,icarbon), &
                     circ_class_n(ipts,ivm,:), ivm, pipe_tune2(ipts,ivm))

             qm_dia_up_half = wood_to_qmdia_up_half(&
                     circ_class_biomass(ipts,ivm,:,:,icarbon), &
                     circ_class_n(ipts,ivm,:), ivm, pipe_tune2(ipts,ivm))


             ! The following is specific for each forest management strategy
             IF(forest_managed(ipts,ivm).EQ.ifm_thin .OR. &
                forest_managed(ipts,ivm).EQ.ifm_uneven ) THEN

                ifm = forest_managed(ipts,ivm)
                ! Forests can only be thinned if it is profitable for the 
                ! managers to do so. This profitability is based on the 
                ! average height of the 100 tallest trees per hectare.  
                ! If it's profitable, the stand is thinned to a target 
                ! relative density index (RDI), which is the ratio of the 
                ! true density to the self thinning density. This target 
                ! depends on the management type, but it is always less 
                ! than unity. Table 1 in Bellassen et al, Ecological 
                ! Modeling (2010) gives values for some of these parameters.

                ! First, calculate the average height of the 100 tallest 
                ! trees on the hectare. Remember that all our trees are the 
                ! same height in a given circ class. Since we don't explicitly 
                ! track trees, we just start at the largest circ classs and 
                ! add trees to our average until we have 100 (Note that this 
                ! number has been externalised as ::ntress_profit)
                height(ipts,ivm,:)=&
                     wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                     ivm,pipe_tune2(ipts,ivm))

                ! Debug
                IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                   ivm == test_pft)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, check for thinning', ivm
                   ! just for printing out some information
                   dia(ipts,ivm,:)= &
                        wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,&
                        pipe_tune2(ipts,ivm))
                   DO icir=1,ncirc     
                      WRITE(numout,*) circ_class_n(ipts,ivm,icir)*m2_to_ha,&
                           height(ipts,ivm,icir), dia(ipts,ivm,icir)
                   ENDDO
                   CALL flush(numout)
                ENDIF

                !-
                ! condition introduced in order to cut the stand at the end
                ! of a spinup run. Put forced_clear_cut=n to avoid
                ! any forced clear-cut 
                IF( forced_clear_cut ) THEN
                      WRITE(numout,*) 'A clear-cut at the begenning of &
                                          the simulation is scheduled'
                         CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                              circ_class_kill(ipts,ivm,:,ifm,icut_clear))
                     CYCLE pft
                ENDIF

                ! Select the 100 tallest treesb and calculate their 
                ! average height
                total_trees=zero
                ave_tree_height=zero
                tree:DO icir=ncirc,1,-1
                   trees=circ_class_n(ipts,ivm,icir)*m2_to_ha
                   IF(trees .GE. ntrees_profit)THEN
                      ave_tree_height=height(ipts,ivm,icir)
                      total_trees=total_trees+trees
                      EXIT tree
                   ELSEIF((trees + total_trees) .GE. ntrees_profit)THEN
                      ave_tree_height=ave_tree_height*total_trees/&
                           REAL(ntrees_profit,r_std)+&
                           (REAL(ntrees_profit,r_std)-total_trees)/&
                           REAL(ntrees_profit,r_std)*&
                           height(ipts,ivm,icir)
                      total_trees=total_trees+trees
                      EXIT tree
                   ELSE
                      ave_tree_height=ave_tree_height*total_trees/&
                           (total_trees+trees)+&
                           trees/(total_trees+trees)*height(ipts,ivm,icir)
                      total_trees=total_trees+trees
                   ENDIF

                ENDDO tree

                ! in order to thin, we need to know the relative 
                ! density index (RDI) of our stand.  This is the 
                ! ratio of the current stand density to
                ! the self-thinning density, Nmax.
                diameters_temp(:)= &
                     wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                     ivm,pipe_tune2(ipts,ivm))

                ! with these diameters, we can check to see if we 
                ! want to clearcut. If the diameters get too big, 
                ! harvesting equipment can't handle the trees anymore. 
                ! On the other hand, we need to make sure there
                ! are enough big trees that it is worthwhile to cut 
                ! them down.
                total_trees=zero
                DO icir=1,ncirc
                   IF(diameters_temp(icir) .GT. largest_tree_dia(ivm)) THEN
                      total_trees = total_trees + circ_class_n(ipts, ivm, icir) * &
                           m2_to_ha
                   ENDIF
                ENDDO

                tot_trees(ipts,ivm) = total_trees

                ! Debug
                IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                    ivm == test_pft)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, testing diameter criterion', ivm
                   WRITE(numout,*) 'total_trees: ', total_trees
                   WRITE(numout,*) 'dens_target values for FM1, FM2, FM3, FM4, FM5: ', dens_target(ivm, ifm_none), dens_target(ivm, &
                   ifm_thin), dens_target(ivm, ifm_cop), dens_target(ivm, ifm_src), dens_target(ivm, ifm_uneven)
                   WRITE(numout,*) 'largest_tree_dia: ', largest_tree_dia(ivm)
                   WRITE(numout,*) 'diameters_temp: ', diameters_temp(:)
                   WRITE(numout,*) 'circ_class_n (in trees per hectare): ', circ_class_n(ipts,ivm,:)*m2_to_ha
                   CALL flush(numout)
                ENDIF
                !-


                ! CONDITION 1: there are less than [dens_target] trees left

                ! Guillaume M. -- CONDITION 1 is disabled under ok_min_density_reset=n:
                ! density alone cannot tell a stand that is open BY DESIGN from a collapsed
                ! one, and a class 1 established at a low PRESCRIBE_RDI_FRAC sits under
                ! dens_target and was clear-cut within the year.
                !
                ! /!\ The enclosing parenthesis is MANDATORY: .AND. binds tighter than .OR.,
                ! so "gate .AND. A .OR. B .OR. C" would gate A only and silently leave the
                ! thin and uneven branches active -- a half-disabled reset.
                !
                ! Guillaume M. -- The discriminator is SIZE, not count: an open young stand
                ! has few THIN stems, a collapsed one has few THICK stems.
                ! age_class_bound(1) is the model's own young/not-young frontier, so no new
                ! parameter is needed.
                !
                ! Guillaume M. -- That size discriminator fails on a crop-tree stand: few
                ! stems ALREADY THICK is the signature it reads as collapse. Under
                ! ok_min_density_rdi the test becomes relative -- a stand held at its RDI
                ! target is sparse by design. Read from rdi_stand: rdi is computed further
                ! down. Not allocated -> keep the historical behaviour, do not disarm a guard.
                cond1_collapsed = .TRUE.
                IF (ok_min_density_rdi .AND. ALLOCATED(rdi_stand)) THEN
                   cond1_collapsed = (rdi_stand(ipts,ivm) .LT. min_density_rdi_floor)
                ENDIF

                IF( ok_min_density_reset .AND. cond1_collapsed .AND. &
                    (qm_dia .GE. age_class_bound(1,ivm)) .AND. ( &
                    ((fm_change_map(ipts,ivm) .EQ. ifm_none)&
                       .AND. SUM(circ_class_n(ipts,ivm,:))*m2_to_ha .LT. dens_target(ivm, ifm_none)) .OR.&
                     ((fm_change_map(ipts,ivm) .EQ. ifm_thin)&
                       .AND. SUM(circ_class_n(ipts,ivm,:))*m2_to_ha .LT. dens_target(ivm, ifm_thin)) .OR.&
                     ((fm_change_map(ipts,ivm) .EQ. ifm_uneven)&
                       .AND. SUM(circ_class_n(ipts,ivm,:))*m2_to_ha .LT. dens_target(ivm, ifm_uneven)) ) ) &
                        THEN
                
                   whychange(ipts,ivm) = 1
                   
                   IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                           WRITE(numout,*) 'CONDITION 1 : ', ivm, whychange(ipts,ivm), ' (pft, whychange)',&
                            ', FM_change_map : ', fm_change_map(ipts,ivm),&
                            ', CCIND : ',  SUM(circ_class_n(ipts,ivm,:))*m2_to_ha, &
                            ', dens_target FM1 FM2 FM5 : ', dens_target(ivm,ifm_none),&
                            dens_target(ivm,ifm_thin), dens_target(ivm,ifm_uneven)
                           CALL flush(numout)
                   ENDIF

                   ! These are conditions which are suitable for switching
                   ! from a management to a conservation strategy. 
                   ! Such a switch is not internalized in the code and 
                   ! needs to be prescribed. The way to prescribe it is 
                   ! by using species/management changes and prescribing 
                   ! the desired management for some of all pixels
                   ! as unmanaged.
                   IF(ok_change_species)THEN

                      ! The criteria for a clearcut are met but check 
                      ! whether we want to switch to a conservation 
                      ! strategy on this site.
                      IF(fm_change_map(ipts,ivm) .EQ. ifm_none) THEN

                         whychange(ipts,ivm) = 1.5

                         ! Change the management strategy (most likely 
                         ! high stand management) to conservation
                         forest_managed(ipts,ivm) = ifm_none

                         IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                 WRITE(numout,*) 'CONDITION 1bis : ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                 ' (pft, whychange, FM)'
                                 CALL flush(numout)
                         ENDIF

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'species change ',&
                                 '- conservation strategy, ',&
                                 'density becoming too low'
                         ENDIF  
                         !-             

                      ELSE

                         ! The species/management change scenario does
                         ! not prescribe a change to a conservation strategy
                         ! but the criteria for a clearcut are met 
                         CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                              circ_class_kill(ipts,ivm,:,ifm,icut_clear))

                      ENDIF
                                
                   ELSE

                      ! Species/management changes are not considered. The 
                      ! criteria for a clearcut are met 
                      CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                           circ_class_kill(ipts,ivm,:,ifm,icut_clear))

                   ENDIF

                   CYCLE pft

                ENDIF
                   

                ! CONDITION 2: There are more than 10 trees for which
                ! the diameter exceeds the threshold, time to cut.
                !
                ! Guillaume M. -- Cut criterion: DIAMETER ONLY. The age criterion and its
                ! OK_AGE_ROTATION flag were removed as redundant with the progressive-harvest
                ! cadence, which already carries the target rotation (see the pre-pass).
                ! DIA_ROTATION_TOL is KEPT: it places the cut INSIDE the last age class, the
                ! bounds being AGE_CLASS_FRAC*LTD with AGE_CLASS_FRAC(nagec-1) below 1-tol.
                !
                ! Guillaume M. -- /!\ The management-intensity maps are ALLOCATED only under
                ! ok_management_intensity, so the ALLOCATED test is mandatory: the OFF path
                ! would otherwise dereference an unallocated array -> SIGSEGV.
                dfac_loc = un
                IF (ALLOCATED(dia_factor)) dfac_loc = dia_factor(ipts)
                ! Guillaume M. -- /!\ Under progressive harvest the diameter test is REMOVED:
                ! being in the LAST AGE CLASS is enough, so the cut diameter acts in a single
                ! place -- the lower bound of that class, (1-DIA_ROTATION_TOL)*
                ! largest_tree_dia. A cold start sets the class diameters without matching
                ! them to the bounds, which otherwise leaves those slots ineligible for decades.
                !
                ! Guillaume M. -- /!\ Gated on ok_progressive_harvest, and this is NOT
                ! cosmetic: the f_cut cadence is what limits HOW MUCH is cut. Without it the
                ! diameter test is the ONLY brake, and removing it would empty every
                ! last-class slot every year. Default .FALSE. -> other configurations keep
                ! the historical behaviour.
                IF (ok_progressive_harvest) THEN
                   ph_dia_ok = .TRUE.
                ELSE
                   ph_dia_ok = (qm_dia_up_half .GE. &
                        largest_tree_dia(ivm)*dfac_loc*(un-dia_rotation_tol))
                ENDIF
                IF( ph_dia_ok .AND. &
                   ifm .EQ. ifm_thin .AND. allow_mgmt_cut) THEN

                   whychange(ipts,ivm) = 2

                   ! Guillaume M. -- PROGRESSIVE_HARVEST: the TRIGGER above is unchanged, it
                   ! still decides WHICH slots are eligible. Only the ACTION becomes
                   ! fractional: ph_f_cut of the slot area is withdrawn, the remainder keeps
                   ! its age, diameter and biomass and is eligible again next year.
                   ! ph_cut_frac hands off to sapiens_kill and age_class_distr.
                   !
                   ! Guillaume M. -- /!\ Snap to the whole-slot clearcut when almost nothing
                   ! would survive: the split restores the surviving density by dividing by
                   ! (1-f), so f near 1 is ill-conditioned on a sliver of ground and is
                   ! physically indistinguishable from cutting the lot.
                   ph_active = .FALSE.
                   IF (ok_progressive_harvest) THEN
                      IF (ph_f_cut_pft(ipts,ivm) .GT. zero .AND. &
                           (un - ph_f_cut_pft(ipts,ivm)) .GT. ph_min_survivor) THEN
                         ph_active = .TRUE.
                         ph_cut_frac(ipts,ivm) = ph_f_cut_pft(ipts,ivm)
                      ENDIF
                   ENDIF

                   IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                           WRITE(numout,*) 'CONDITION 2 : ', ivm, whychange(ipts,ivm), ' (pft, whychange)',&
                         ', qm_dia_up_half : ', qm_dia_up_half,&
                         ', largest_tree_dia : ', largest_tree_dia(ivm)
                           CALL flush(numout)
                   ENDIF

                   ! Debug
                   IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                      WRITE(numout,*) 'Sapiens_forestry_main, diameter exceeds ',&
                          'treshold, clearcut', ivm
                      WRITE(numout,*) 'total_trees: ', total_trees
                      WRITE(numout,*) 'dens_target FM1, FM2, FM3, FM4, FM5: ', dens_target(ivm, ifm_none),&
                              dens_target(ivm, ifm_thin), dens_target(ivm, ifm_cop),&
                              dens_target(ivm, ifm_src), dens_target(ivm, ifm_uneven)
                      WRITE(numout,*) 'largest_tree_dia: ', &
                           largest_tree_dia(ivm)
                      WRITE(numout,*) 'diameters_temp: ', diameters_temp(:)
                      WRITE(numout,*) 'circ_class_n (trees per hectare): ', &
                           circ_class_n(ipts,ivm,:)*m2_to_ha
                      CALL flush(numout)
                   ENDIF
                   !-

                   ! These are conditions which are suitable for switching
                   ! from a management to a conservation strategy. Such a 
                   ! switch is not internalized in the code and needs to be 
                   ! prescribed. The way to prescribe it is by using 
                   ! species/management changes and prescribing the desired 
                   ! management for some of all pixels as unmanaged.
                   IF(ok_change_species)THEN

                      ! The criteria for a clearcut are met but check whether 
                      ! we want to switch to a conservation strategy on this 
                      ! site.
                      IF( fm_change_map(ipts,ivm) .EQ. ifm_none) THEN

                         whychange(ipts,ivm) = 2.5

                         ! Change the management strategy (most likely 
                         ! high stand management) to conservation
                         forest_managed(ipts,ivm) = ifm_none

                         IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                 WRITE(numout,*) 'CONDITION 2bis : ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                ' (pft, whychange, FM)'
                                 CALL flush(numout)
                         ENDIF

                         ! debug
                         IF(printlev_loc>=4)THEN
                             WRITE(numout,*) 'species change - ',&
                                'conservation strategy, ',&
                                'exceeded largest diameter'
                         ENDIF
                         !-

                      ELSE

                         ! The species/management change scenario does
                         ! not prescribe a change to a conservation strategy
                         ! but the criteria for a clearcut are met
                         IF (ph_active) THEN
                            CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                                circ_class_kill(ipts,ivm,:,ifm,icut_clear), &
                                frac=ph_f_cut_pft(ipts,ivm))
                         ELSE
                            CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                                circ_class_kill(ipts,ivm,:,ifm,icut_clear))
                         ENDIF

                      ENDIF

                   ELSE

                      ! Species/management changes are not considered. The
                      ! criteria for a clearcut are met
                      IF (ph_active) THEN
                         CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                             circ_class_kill(ipts,ivm,:,ifm,icut_clear), &
                             frac=ph_f_cut_pft(ipts,ivm))
                      ELSE
                         CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                             circ_class_kill(ipts,ivm,:,ifm,icut_clear))
                      ENDIF

                   ENDIF

                   CYCLE pft
                ENDIF

                ! CONTION 3: a clearcut is expected if the trees start
                ! growing too slowly. This is measured by looking at 
                ! the mean annual increment, which is simply the trunk 
                ! volume divided by the age of the tree.  We compare 
                ! this to average tree growth over the past five years 
                ! (i.e. the tree volume divided by 5), called the 
                ! periodic annual increment.  If the PAI is smaller than 
                ! the MAI, the trees are past their prime and the owner 
                ! will cut them and plant anew.

                ! We calculate the MAI and the PAI at the same place where 
                ! we increment the age of the stand, so here we just have 
                ! to compare them.  We need to check that we have enough 
                ! years from the last cut to calculate a good PAI.
                IF(PAI(ipts,ivm) .LT. MAI(ipts,ivm) .AND. &
                     last_cut(ipts,ivm) .GT. n_pai .AND. allow_mgmt_cut) THEN

                    whychange(ipts,ivm) = 3

                    IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                            WRITE(numout,*) 'CONDITION 3 : ', ivm, whychange(ipts,ivm), ' (pft, whychange)',&
                            ', PAI, MAI, last_cut, n_pai : ', PAI(ipts,ivm), MAI(ipts,ivm), last_cut(ipts,ivm), n_pai
                            CALL flush(numout)
                    ENDIF

                    ! Debug
                    IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                        ivm == test_pft)THEN
                       WRITE(numout,*) 'Sapiens_forestry_main, harvesting ',&
                            'increment PAI<MAI'
                       WRITE(numout,*) 'ipts, ivm, mai, pai: ',ipts, &
                            ivm, mai(ipts,ivm), pai(ipts,ivm)
                    ENDIF
                    !-

                    IF(ok_change_species)THEN

                       ! The criteria for a clearcut are met but check 
                       ! whether we want to switch to a conservation 
                       ! strategy on this site.
                       IF( fm_change_map(ipts,ivm) .EQ. ifm_none) THEN

                          whychange(ipts,ivm) = 3.5

                          ! Change the management strategy (could be 
                          ! high stand, coppice and short rotation 
                          ! coppice) to conservation
                          forest_managed(ipts,ivm) = ifm_none

                          IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                  WRITE(numout,*) 'CONDITION 3bis : ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                ' (pft, whychange, FM)'
                                  CALL flush(numout)
                          ENDIF

                          IF(printlev_loc>=4)THEN
                             WRITE(numout,*) 'species change - ',&
                                  'conservation strategy, ', &
                                  'stand growth is too low'
                          ENDIF  

                      ELSE

                         ! The species/management change scenario does
                         ! not prescribe a change to a conservation strategy
                         ! but the criteria for a clearcut are met 
                         CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                              circ_class_kill(ipts,ivm,:,ifm,icut_clear))

                      ENDIF
                                 
                   ELSE

                      ! Species/management changes are not considered. The 
                      ! criteria for a clearcut are met
                      CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                           circ_class_kill(ipts,ivm,:,ifm,icut_clear))
                        
                   ENDIF

                  CYCLE pft
                ENDIF

                ! The self-thinning and yield relationships are for 
                ! diameters expressed in cm
                rdi(ipts,ivm) = (SUM(circ_class_n(ipts,ivm,:))*m2_to_ha)/ &
                        Nmax(qm_dia*m_to_cm, alpha_self_thinning(ipts,ivm), ivm, ifm) 

                !CONDITION 4
                IF(fm_change_map(ipts,ivm)==5 .AND. &
                   rdi(ipts,ivm) <= rdi_max(ivm,ifm_uneven)*1.1 .AND. &
                   rdi(ipts,ivm) >= rdi_max(ivm,ifm_uneven)*0.9) THEN
                  forest_managed(ipts,ivm) = fm_change_map(ipts,ivm)

                  whychange(ipts,ivm) = 4

                  IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                          WRITE(numout,*) 'CONDITION 4 : ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                          ' (pft, whychange, FM)',&
                          ', rdi, rdi_max : ', rdi(ipts,ivm), rdi_max(ivm,ifm_uneven)
                          CALL flush(numout)
                  ENDIF

                ENDIF

                !CONDITION 5

                IF(fm_change_map(ipts,ivm)==1 .AND. &
                   rdi(ipts,ivm) >= rdi_max(ivm,ifm_none)*0.9) THEN
                  forest_managed(ipts,ivm) = fm_change_map(ipts,ivm)

                  whychange(ipts,ivm) = 5

                  IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                          WRITE(numout,*) 'CONDITION 5 : ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                          ' (pft, whychange, FM)',&
                          ', rdi, rdi_max : ', rdi(ipts,ivm), rdi_max(ivm,ifm_none)
                          CALL flush(numout)
                  ENDIF

                ENDIF

                
                ! rdi_target_upper (based on the yield tables)
                ! rdi depends on the observed self-thinning 
                ! relationship and rdi_target_upper
                ! Polynomial approach
                ifm=forest_managed(ipts,ivm)
                CALL calculate_rdi_boundaries(qm_dia_up_half, ivm, ifm, &
                                             rdi_target_upper(ipts,ivm), &
                                             rdi_target_lower(ipts,ivm))

                ! Guillaume M. -- /!\ POLYNOMIAL SHORT-CIRCUIT for the INTERMEDIATE classes.
                ! The polynomials give the target RDI from the DIAMETER, tuned per PFT and
                ! per management TYPE but blind to the intensity gradient: at a given
                ! diameter the curve caps below the observed stocking and the stand can
                ! never densify.
                !
                ! Guillaume M. -- Functional split of the classes:
                !   class 1        : RDI FREE  -> post-disturbance opening, edge signal;
                !   classes 2..n-1 : RDI FIXED -> growth phase, management acts HERE;
                !   class nagec    : RDI FREE  -> harvest phase, overstocking regulated by
                !                    the disturbances (thinning already suppressed there).
                !
                ! Guillaume M. -- OK_RDI_BAND_CLASS1 moves class 1 into the banded interval.
                ! Leaving it RDI-free traps a dense stand: it can neither be thinned nor grow
                ! past age_class_bound(1), so it never graduates. Same bound is applied to
                ! thin_this_slot below -- a band never read, or a thinning against the
                ! polynomial, would each be worse than the trap.
                !
                ! Guillaume M. -- RDI being RELATIVE (N/Nmax(qmd, alpha)), a fixed value does
                ! not erase species differences: they pass through alpha in the denominator.
                ! Same gate as the suppression of thinning in the last class: both express
                ! the FUNCTIONAL-class paradigm and must not be activated separately.
                iage_thin_lo = 2
                IF (ok_rdi_band_class1) iage_thin_lo = 1
                IF (ok_clearcut_last_class .AND. nagec > 1) THEN
                   iage_th = ivm - start_index(agec_group(ivm)) + 1
                   IF (iage_th >= iage_thin_lo .AND. iage_th <= nagec_pft(agec_group(ivm)) - 1) THEN
                      rdi_target_upper(ipts,ivm) = rdi_band_high(ivm)
                      rdi_target_lower(ipts,ivm) = rdi_band_low(ivm)

                      ! Guillaume M. -- Density is in TRANSITION through the improvement
                      ! class, where the crop trees are singled out. A constant target set at
                      ! its entrance is a STEP, not a thinning -- a clearcut under another
                      ! name. Progress is read on the DIAMETER, which already defines class
                      ! membership and carries the pixel's fertility.
                      !
                      ! Guillaume M. -- GEOMETRIC interpolation: a ratio is interpolated and
                      ! stem number decays multiplicatively -- a constant removal RATE, not a
                      ! constant subtraction. The band keeps its relative width so the trigger
                      ! follows the target down instead of drifting away from it.
                      IF (ok_rdi_ramp .AND. &
                           iage_th .EQ. nagec_pft(agec_group(ivm)) - 1) THEN
                         ramp_u = (qm_dia - age_class_bound(iage_th-1,ivm)) &
                              / MAX(age_class_bound(iage_th,ivm) &
                                    - age_class_bound(iage_th-1,ivm), min_stomate)
                         ramp_u = MIN(un, MAX(zero, ramp_u))
                         ramp_lo = rdi_band_low(ivm)**(un - ramp_u) &
                              * rdi_ramp_out(ivm)**ramp_u
                         rdi_target_lower(ipts,ivm) = ramp_lo
                         rdi_target_upper(ipts,ivm) = ramp_lo &
                              * rdi_band_high(ivm) / MAX(rdi_band_low(ivm), min_stomate)
                      ENDIF
                   ENDIF
                ENDIF

                ! Guillaume M. -- Publish the RDI target management ACTUALLY uses. The field
                ! RDI_TARGET_LOWER cannot serve: stomate_lpj zeroes it and recomputes it from
                ! the polynomial WITHOUT the band override above. Recorded here rather than at
                ! the CALL -- nothing writes rdi_target_lower in between, so this is the value
                ! handed to thinning, and an unthinned slot still shows its target.
                IF (.NOT. ALLOCATED(rdi_thin_limit_diag)) THEN
                   ALLOCATE(rdi_thin_limit_diag(npts,nvm))
                   rdi_thin_limit_diag(:,:) = undef
                ENDIF
                rdi_thin_limit_diag(ipts,ivm) = rdi_target_lower(ipts,ivm)

                ! Debug
                IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, PFT ',ivm
                   WRITE(numout,*) 'rdi', rdi(ipts,ivm)
                   WRITE(numout,*) 'rdi_target_upper calcul', &
                        rdi_target_upper(ipts,ivm)
                   WRITE(numout,*) 'dom_dia', qm_dia_up_half*m_to_cm 
                   WRITE(numout,*) 'dia', qm_dia*m_to_cm 
                   WRITE(numout,*) 'Sapiens_forestry_main, rdi print ',rdi(ipts,ivm),&
                        rdi_target_upper(ipts,ivm),&
                        rdi_target_lower(ipts,ivm) 
                   WRITE(numout,*) "Nmax ",&
                           Nmax(qm_dia*m_to_cm, alpha_self_thinning(ipts,ivm), ivm, ifm)
                ENDIF
                !-

                ! we only need to thin if the rdi of our stand is outside
                ! the acceptable rdi range given by rdi_target_upper and
                ! rdi_target_lower
                !
                ! Guillaume M. -- /!\ NO THINNING IN THE LAST AGE CLASS. Classes 1..nagec-1
                ! are the GROWTH phase, where thinning is the silvicultural tool and where
                ! management intensity acts; the last class is the HARVEST phase -- it is
                ! cut, not thinned. Thinning therefore happens ONLY in the intermediate
                ! classes, exactly where the fixed RDI band above is defined.
                !
                ! Guillaume M. -- The mechanism forces it: beyond r_max*LTD thinning switches
                ! to from-above (-thinstrat) and removes the LARGEST trees. The last class
                ! spans [1-DIA_ROTATION_TOL ; 1]*LTD, hence is ALWAYS above that threshold,
                ! so every mature stand had its upper tail cropped -- a direct brake on
                ! diameter in the one class where it must grow.
                !
                ! Guillaume M. -- /!\ Without thinning AND without an RDI cap here, only
                ! self-thinning and the EXPLICIT disturbances regulate density: the external
                ! acceptance criterion is standing biomass against inventories. Gated on
                ! ok_clearcut_last_class, which already confines the management cut to this
                ! class -- without it the class would lose both its outlet and its regulation.
                thin_this_slot = .TRUE.
                IF (ok_clearcut_last_class .AND. nagec > 1) THEN
                   iage_th = ivm - start_index(agec_group(ivm)) + 1
                   thin_this_slot = (iage_th >= iage_thin_lo .AND. &
                        iage_th <= nagec_pft(agec_group(ivm)) - 1)
                ENDIF
                ! Guillaume M. -- FORCED_THINNING: on the young maturity classes the thinning
                ! no longer waits for the RDI band to be crossed, it happens EVERY year at a
                ! fixed intensity. The band is calibrated on the temperate self-thinning line
                ! while ref_alpha_self_thin is calibrated per MTC: where the line tolerates
                ! more stems the RDI stays under the band and the stand is never thinned.
                forced_thin_slot = .FALSE.
                IF (ok_forced_thinning .AND. nagec > 1) THEN
                   iage_ft = ivm - start_index(agec_group(ivm)) + 1
                   forced_thin_slot = (iage_ft >= 1 .AND. iage_ft <= forced_thin_last_class)
                ENDIF

                ! Guillaume M. -- The forced rule REPLACES the RDI trigger on those slots
                ! rather than adding to it: applying both in the same year would compound two
                ! thinnings. The target is expressed in RDI so that `thinning` below is reused
                ! unchanged, keeping the from-below/from-above switch and the icut_thin path.
                rdi_thin_goal = rdi_target_lower(ipts,ivm)
                IF (forced_thin_slot .AND. rdi(ipts,ivm) .GT. min_stomate) THEN
                   rdi_thin_goal = rdi(ipts,ivm) * (un - forced_thin_intensity(ivm))
                ENDIF

                IF((forced_thin_slot .AND. rdi(ipts,ivm) .GT. min_stomate) .OR. &
                     (.NOT. forced_thin_slot .AND. &
                      rdi(ipts,ivm) .GT. rdi_target_upper(ipts,ivm) .AND. thin_this_slot)) THEN

                   ! we need the circumference of each model tree in meters 
                   circ_temp(:)= wood_to_circ(&
                        circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,&
                        pipe_tune2(ipts,ivm))

                   ! Debug
                   IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                        ivm == test_pft)THEN
                      DO icir=1,ncirc
                         WRITE(numout,*) 'Sapiens_forestry_main thinning now: ',&
                              ipts,ivm,icir
                         WRITE(numout,*) 'Sapiens_forestry_main height, ccn, ',&
                              height(ipts,ivm,icir),&
                              circ_class_n(ipts,ivm,icir)*m2_to_ha
                      ENDDO
                      WRITE(numout,*) "circ_temp ",circ_temp(:)
                      WRITE(numout,*) "NOTE: assumes ncirc = 3"
                      WRITE(numout,*) "circ_class_biomass(1) ",&
                           SUM(circ_class_biomass(ipts,ivm,1,:,icarbon))
                      WRITE(numout,*) "circ_class_biomass(2) ",&
                           SUM(circ_class_biomass(ipts,ivm,2,:,icarbon))
                      WRITE(numout,*) "circ_class_biomass(3) ",&
                           SUM(circ_class_biomass(ipts,ivm,3,:,icarbon))
                   ENDIF
                   !-

                   ! As long as the trees are relatively small (=66% of their
                   ! target diameter), we thin from below. If the trees are
                   ! getting larger, thinning will be from above. This switch
                   ! in thinning regime is essential to maintain a realistic
                   ! BA for mature/older stands.
                   if (qm_dia_up_half >= (largest_tree_dia(ivm)*r_max)) then
                     CALL thinning(taumin(ivm),taumax(ivm),-(thinstrat(ivm)),&
                        diameters_temp(:),circ_temp(:),&
                        circ_class_kill(ipts,ivm,:,ifm,icut_thin),&
                        circ_class_n(ipts,ivm,:),rdi_thin_goal,&
                        ivm,ipts, ifm) 
                   else
                     CALL thinning(taumin(ivm),taumax(ivm),thinstrat(ivm),&
                        diameters_temp(:),circ_temp(:),&
                        circ_class_kill(ipts,ivm,:,ifm,icut_thin),&
                        circ_class_n(ipts,ivm,:),rdi_thin_goal,&
                        ivm,ipts, ifm)
                   endif

                   ! Debug
                   IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                       ivm == test_pft)THEN
                      WRITE(numout,*) 'Sapiens_forestry_main, circ_class_kill ',&
                           circ_class_kill(ipts,ivm,:,ifm,icut_thin)
                      WRITE(numout,*) 'circ_class_n(ipts,ivm,:) ',&
                           circ_class_n(ipts,ivm,:)
                      DO icir=1,ncirc
                         WRITE(numout,*) 'Sapiens_forestry_main, circ_biomass(iparts)',&
                              icir,circ_class_biomass(ipts,ivm,icir,:,icarbon)
                      ENDDO
                   ENDIF
                   !-

                ENDIF ! rdi .GT. rdi_target_upper

             ! Coppicing
             ELSEIF(forest_managed(ipts,ivm) .EQ. ifm_cop) THEN

                ! This is for coppicing.
                ! It makes no sense to coppice an evergreen tree species.
                IF(is_needleleaf(ivm)) THEN

                   WRITE(numout,*) 'WARNING', is_needleleaf(ivm)
                   WRITE(numout,*) 'Sapiens_forestry_main, Coppicing (FM = 3) should  ' &
                        // 'only be done on deciduous trees.'
                   WRITE(numout,*) 'ipts,ivm,forest_managed(ipts,ivm): ',&
                        ipts,ivm,forest_managed(ipts,ivm) 
                   IF (err_act.GT.1) CALL ipslerr_p(2,'sapiens_forestry_main', &
                        'Don t coppice conifers','','')
                   
                ENDIF

                ! Depending on the age of the stand we have different options.
                ! The current system is to do nothing until the diameter of the
                ! trees is greater than a certain value.  Once this happens, we
                ! coppice, saving the number of trees that exist in this stand.
                ! The trees then regrow until they reach this diameter again.
                ! Some trees will have died due to mortality.  At the second
                ! (and subsequent) cuts, we redistribute the below ground root 
                ! mass to have the same number of trees as the first cut, 
                ! killing the fine roots and harvesting aboveground biomass.

                ! Check first to see if the average diameter of the
                ! stand is greater than the prescribed treshold parameter.
                diameters_temp(:)=wood_to_dia(&
                     circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,&
                     pipe_tune2(ipts,ivm))

                rdi(ipts,ivm) = (SUM(circ_class_n(ipts,ivm,:))*m2_to_ha)/ &
                        Nmax(qm_dia*m_to_cm, alpha_self_thinning(ipts,ivm), ivm, ifm)

                
                ! Debug
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, Do we coppice? '
                   WRITE(numout,*) 'ipts,ivm,SUM(circ_class_n(ipts,ivm,:)) ',&
                        ipts,ivm,SUM(circ_class_n(ipts,ivm,:))
                ENDIF
                !-

                ! We have already coppiced once, so we know our target
                ! tree density.
                IF(qm_dia_up_half  .GE. coppice_diameter(ivm) .AND. &
                     coppice_dens(ipts,ivm) .GT. zero)THEN

                   ! Debug
                   IF(printlev_loc>=4) WRITE(numout,*) 'Sapiens_forestry_main, '&
                        // 'Coppicing second cut. '
                   !-

                   IF(ok_change_species)THEN

                      ! The criteria coppicing are met but check whether 
                      ! we want to switch to a conservation strategy on 
                      ! this site.
                      IF(fm_change_map(ipts,ivm) .EQ. ifm_none) THEN

                         ! Change the management strategy in this case 
                         ! coppicing to conservation
                         forest_managed(ipts,ivm) = ifm_none

                         whychange(ipts,ivm) = 6

                         IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                 WRITE(numout,*) 'CONDITION 6: ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                 ' (pft, whychange, FM)',&
                                 ', qm_dia_up_half, coppice_diameter, coppice_dens : ',&
                                 qm_dia_up_half, coppice_diameter(ivm), coppice_dens(ipts,ivm)
                                 CALL flush(numout)
                         ENDIF


                      ELSEIF((fm_change_map(ipts,ivm).EQ.ifm_thin) .OR. &
                         (fm_change_map(ipts,ivm).EQ.ifm_uneven) .OR. & 
                         (fm_change_map(ipts,ivm).EQ.ifm_src)) THEN    

                         ifm = fm_change_map(ipts,ivm) 
                         ! We would like to change the forest management so 
                         ! this is a good opportunity. However, we can not 
                         ! simply coppice because then the roots will stay
                         ! on the site. We have two options depending on
                         ! whether we also want to change the species or
                         ! not.
                         
                         ! Calculate the pft corresponding with the youngest
                         ! age class of this species.
                         igroup_new = species_change_map(ipts,ivm)
                         igroup_old = agec_group(ivm)

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'replant coppice 1+?'
                            WRITE(numout,*) 'igroup_new, igroup_old, ivm, ', &
                                 igroup_new, igroup_old, ivm
                         ENDIF
                         !-
                         
                         IF(igroup_new.EQ.igroup_old)THEN

                            ! We have the intention to replant with the same
                            ! species and change to high stand management so
                            ! it is better not to harvest the trees and let
                            ! the coppice develop in standards.
                            forest_managed(ipts,ivm) = ifm

                            whychange(ipts,ivm) = 7

                            IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                    WRITE(numout,*) 'CONDITION 7: ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                    ' (pft, whychange, FM)',&
                                    ', igroup_new, igroup_old : ',&
                                    igroup_new, igroup_old
                                    CALL flush(numout)
                            ENDIF

                         ELSE

                            ! We want to replace this coppice by another 
                            ! species and to change the management strategy.
                            ! We have to make sure that the whole stand is
                            ! harvested. If so biomass will be zero and
                            ! the sapiens_forestry_flag_species_change routine will flag this
                            ! PFT to be replanted in species_change the
                            ! management will be changed
                            CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                              circ_class_kill(ipts,ivm,:,ifm,icut_clear))

                            whychange(ipts,ivm)=8
                            IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                    WRITE(numout,*) 'CONDITION 8: ', ivm, whychange(ipts,ivm), ' (pft, whychange)'
                                    CALL flush(numout)
                            ENDIF

                         ENDIF !species_change

                      ELSE

                         ! The species/management change scenario does
                         ! not prescribe a change to a conservation strategy
                         ! but the criteria for coppicing are met. We have 
                         ! already coppiced once, so we know our target
                         ! tree density.  We just have to mark all the trees
                         ! for coppicing.
                         circ_class_kill(ipts,ivm,:,ifm_cop,icut_cop2) = & 
                              circ_class_n(ipts,ivm,:)

                      ENDIF
                                 
                   ELSE

                      ! Species/management changes are not considered. The 
                      ! criteria for a coppice are met.  We have already 
                      ! coppiced once, so we know our target tree density.  
                      ! We just have to mark all the trees for coppicing.
                      circ_class_kill(ipts,ivm,:,ifm_cop,icut_cop2) = &
                           circ_class_n(ipts,ivm,:)

                   ENDIF
                   
                  
                ! First time coppicing still need to save the coppice density
                ELSEIF(qm_dia_up_half  .GE. coppice_diameter(ivm))THEN
                   
                   ! debug
                   IF(printlev_loc>=4)THEN
                      WRITE(numout,*) 'Sapiens_forestry_main, Coppicing first cut. '
                      WRITE(numout,*) 'coppice_dens(ipts,ivm): ',&
                           SUM(circ_class_n(ipts,ivm,:))
                   ENDIF

                   ! save this tree density.  We will reallocate to this in 
                   ! future cuts.
                   coppice_dens(ipts,ivm)=SUM(circ_class_n(ipts,ivm,:))

                   IF(ok_change_species)THEN

                      ! The criteria for coppicing are met but check whether 
                      ! we want to switch to a conservation strategy on this 
                      ! site.
                      IF( fm_change_map(ipts,ivm) .EQ. ifm_none) THEN

                         ! Change the management strategy in this case
                         ! coppicing to conservation
                         forest_managed(ipts,ivm) = ifm_none

                         whychange(ipts,ivm) = 9

                         IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                 WRITE(numout,*) 'CONDITION 9: ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                 ' (pft, whychange, FM)',&
                                 ', qm_dia_up_half, coppice_diameter : ',&
                                 qm_dia_up_half, coppice_diameter(ivm)
                                 CALL flush(numout)
                         ENDIF


                      ELSEIF((fm_change_map(ipts,ivm) .EQ. ifm_thin) .OR. &
                         (fm_change_map(ipts,ivm) .EQ. ifm_uneven) .OR. &
                         (fm_change_map(ipts,ivm) .EQ. ifm_src)) THEN

                         ifm = fm_change_map(ipts,ivm)
                         ! We would like to change the forest management so 
                         ! this is a good opportunity. However, we can not 
                         ! simply coppice because then the roots will stay
                         ! on the site. We have two options depending on
                         ! whether we also want to change the species or
                         ! not.
                         
                         ! Calculate the pft corresponding with the youngest
                         ! age class of this species.
                         igroup_new = species_change_map(ipts,ivm)
                         igroup_old = agec_group(ivm)

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'replant coppice 1?'
                            WRITE(numout,*) 'igroup_new, igroup_old, ivm, ', &
                                 igroup_new, igroup_old, ivm
                         ENDIF
                         !-
                         
                         IF(igroup_new.EQ.igroup_old)THEN

                            ! We have the intention to replant with the same
                            ! species and change to high stand management so
                            ! it is better not to harvest the trees and let
                            ! the coppice develop in standards. Stop coppicing
                            ! and change the management
                            forest_managed(ipts,ivm) = ifm

                            whychange(ipts,ivm) = 10

                            IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                    WRITE(numout,*) 'CONDITION 10: ', ivm, whychange(ipts,ivm), forest_managed(ipts,ivm),&
                                    ' (pft, whychange, FM)',&
                                    ', igroup_new, igroup_old : ',&
                                    igroup_new, igroup_old
                                    CALL flush(numout)
                            ENDIF

                         ELSE

                            ! We want to replace this coppice by another 
                            ! species and by another management strategy.
                            ! We have to make sure that the whole stand is
                            ! harvested. If so, biomass will be zero and
                            ! the sapiens_forestry_flag_species_change routine will flag this
                            ! PFT to be replanted. In species_change the
                            ! management will be changed and the new PFT
                            ! will be set.
                            CALL clearcut_harvest(circ_class_n(ipts,ivm,:), &
                              circ_class_kill(ipts,ivm,:,ifm,icut_clear))

                            whychange(ipts,ivm)=11

                            IF (printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft) THEN
                                    WRITE(numout,*) 'CONDITION 11: ', ivm, whychange(ipts,ivm), ' (pft, whychange)'
                                    CALL flush(numout)
                            ENDIF

                         ENDIF !species_change

                      ELSE

                         ! The species/management change scenario does
                         ! not prescribe a change to a conservation strategy
                         ! but the criteria for a coppicing are met. This is
                         ! the first coppice. Mark all the trees
                         ! for coppicing.
                         circ_class_kill(ipts,ivm,:,ifm_cop,icut_cop1) = &
                              circ_class_n(ipts,ivm,:)

                      ENDIF
                      
                   ELSE

                      ! Species/management changes are not considered. The 
                      ! criteria for a coppice are met. Mark all the trees 
                      ! for coppicing.
                      circ_class_kill(ipts,ivm,:,ifm_cop,icut_cop1) = &
                           circ_class_n(ipts,ivm,:)
                      
                   ENDIF

                ELSE IF (.NOT. ok_coppice_2class .AND. &
                     rdi(ipts,ivm) .GT. rdi_max(ivm,ifm_cop)) THEN
                   ! Guillaume M. -- COPPICE_2CLASS: no thinning in a coppice. Reineke's RDI
                   ! describes seed-origin stems competing symmetrically for space; a coppice
                   ! carries several shoots per stool, their number set by the inherited root
                   ! system. The index measures nothing here and rdi_max_cop is calibrated on
                   ! high forest, so this branch compared a meaningless number to a threshold.
                   !
                  ! we need the circumference of each model tree in meters
                   circ_temp(:)= wood_to_circ(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                        ivm,pipe_tune2(ipts,ivm))

                   CALL thinning(taumin(ivm),taumax(ivm),thinstrat(ivm),&
                        diameters_temp(:),circ_temp(:),&
                        circ_class_kill(ipts,ivm,:,ifm,icut_thin),&
                        circ_class_n(ipts,ivm,:),rdi_max(ivm,ifm_cop)*0.9,&
                        ivm,ipts, ifm)
  
                ELSE
                   
                   ! We don't need to coppice anything.

                ENDIF

             ! This is for short rotation coppices.
             ELSEIF(forest_managed(ipts,ivm).EQ.ifm_src)THEN
                
                ! Guillaume M. -- Guard WIDENED. "SRC is limited to poplar and willow" holds
                ! for the North European short-rotation coppice but not for Iberian
                ! eucalyptus, an EVERGREEN broadleaf whose whole silviculture rests on stump
                ! resprouting. FM_SRC_PFT is an explicit user assertion: honoured, but traced
                ! at printlev 2 so the case stays visible.
                !
                ! Guillaume M. -- The coppice machinery does NOT depend on deciduousness: the
                ! cut is age/area accounting on MOD(age_stand, src_rot_length), and
                ! first/second_coppice_cut reference neither phenology, reserves nor
                ! senescence. is_deciduous is only a derivative of pheno_model.
                IF(.NOT. is_deciduous(ivm) .AND. .NOT. fm_src_pft(ivm))THEN

                   WRITE(numout,*) 'ERROR: Sapiens_forestry_main, Short rotation ' &
                        // ' coppicing (FM = 4)' &
                        // 'should only be done on deciduous trees.'
                   WRITE(numout,*) 'ipts,ivm,forest_managed(ipts,ivm): ',&
                       ipts,ivm,forest_managed(ipts,ivm)
                   IF (err_act.GT.1) CALL ipslerr_p(3,'ERROR: sapiens_forestry_main', &
                        'SRC should only be done on deciduous trees.',&
                        'Set FM_SRC_PFT for a resprouting evergreen (e.g. eucalyptus).','')

                ELSEIF(.NOT. is_deciduous(ivm))THEN

                   IF (printlev_loc >= 2) WRITE(numout,*) &
                        'sapiens_forestry_main : taillis sur un SEMPERVIRENT, autorise ' &
                        // 'par FM_SRC_PFT. PFT ', ivm

                ENDIF

                ! Debug
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, Do we SRC? '
                   WRITE(numout,*) 'ipts,ivm,SUM(circ_class_n(ipts,ivm,:)) ',&
                        ipts,ivm,SUM(circ_class_n(ipts,ivm,:))
                   WRITE(numout,*) 'age_stand(ipts,ivm),src_rot_length(ivm),' &
                        // ' coppice_dens(ipts,ivm): ',&
                        age_stand(ipts,ivm),src_rot_length(ivm), &
                        coppice_dens(ipts,ivm)
                   WRITE(numout,*) 'src_nrots(ivm) ',&
                        src_nrots(ivm)
                ENDIF
                !-

                ! NOTE: SRC is not implemented as a historic forest 
                ! management strategy but only used as a future strategy. 
                ! Hence, for the moment we never have the situation where we 
                ! decide we want to convert an SRC into an unmanaged forest. 
                ! This is not a very realist conversion anyway. Contrary to 
                ! ifm_none, ifm_cop and ifm_thin we don not check whether 
                ! the future management strategy is ifm_none.

                ! The qualifications for doing a short rotation coppice are 
                ! fairly simple.  It's based on the age of the stand.  
                ! We harvest all the aboveground biomass.  After a certain 
                ! number of rotations, we harvest the aboveground and kill 
                ! the belowground.
                ! http://www.coppiceresources.co.uk/SRC.asp says that 30 
                ! years between total harvests.
                IF(MOD(age_stand(ipts,ivm),src_rot_length(ivm)) == 0)THEN

                   nrotations=age_stand(ipts,ivm)/src_rot_length(ivm)

                   IF(nrotations == src_nrots(ivm))THEN

                      ! Debug
                      IF(printlev_loc>=4) WRITE(numout,*) 'Sapiens_forestry_main, ' &
                           // 'SRC final harvest. '
                      !-

                      ! This is the final harvest, where we will kill
                      ! all belowground biomass.
                      circ_class_kill(ipts,ivm,:,ifm_src,icut_cop3)=&
                           circ_class_n(ipts,ivm,:)

                   ELSEIF(nrotations == 1)THEN

                      ! Debug
                       IF(printlev_loc>=4) WRITE(numout,*) 'Sapiens_forestry_main, '&
                           // 'SRC first harvest. '
                       !-

                      ! This is the first harvest, which will lead
                      ! to coppices. Save this tree density. We will 
                      ! reallocate to this in future cuts.
                      coppice_dens(ipts,ivm)=SUM(circ_class_n(ipts,ivm,:))
                      
                      ! Mark all trees for the first coppice cut.
                      circ_class_kill(ipts,ivm,:,ifm_src,icut_cop1)= &
                           circ_class_n(ipts,ivm,:)

                   ELSEIF(nrotations < src_nrots(ivm))THEN

                      ! Debug
                      IF(printlev_loc>=4) WRITE(numout,*) 'Sapiens_forestry_main, ' &
                           // 'SRC standard harvest. '
                      !-

                      ! This is a standard harvest where we take the
                      ! aboveground biomass and leave the belowground.
                      ! Mark all trees for a standard coppice cut.
                      circ_class_kill(ipts,ivm,:,ifm_src,icut_cop2)=&
                           circ_class_n(ipts,ivm,:)

                  ELSE

                      WRITE(numout,*) 'ERROR: Sapiens_forestry_main, How did we arrive '&
                           // 'here? We must have '
                      WRITE(numout,*) 'somehow skipped the harvest in SRC. '
                      WRITE(numout,*) 'ipts, ivm: ',ipts,ivm
                      WRITE(numout,*) 'nrotations,src_nrots(ivm): ',&
                           nrotations,src_nrots(ivm)
                      IF(err_act.GT.1) CALL ipslerr_p(3,'ERROR: sapiens_forestry_main', &
                           'Somehow skipped harvest in SRC.',&
                           'Look in the output file for ERROR.','')

                   ENDIF ! check to see how many rotations we've done


                ENDIF ! checking if it's a harvest year

             ELSE

                WRITE(numout,*) 'Sapiens_forestry_main, We do not have any action '&
                     //'for this management strategy'
                WRITE(numout,*) 'forest_managed(ipts,ivm),ipts,ivm ', &
                     forest_managed(ipts,ivm),ipts,ivm
                IF(err_act.GT.1) CALL ipslerr_p(3,'ERROR: sapiens_forestry_main', &
                     'No action for this management strategy.',&
                     'Look in the output file for ERROR.','')

             ENDIF ! test of management strategy

          END IF
          
       ENDDO pft

    ENDDO

    ! Debug
    IF(printlev_loc>=4)THEN
       DO ipts=1,npts
          DO ivm=1,nvm

             DO icir = 1,ncirc
                IF(ipts == test_grid .AND. ivm == test_pft)THEN
                   WRITE(numout,*) 'Sapiens_forestry_main, How many trees do we kill?'
                   WRITE(numout,*) 'Name  icir  test_pft  circ_class_kill  '
                   WRITE(numout,*) 'Thin ',icir,test_pft,&
                     circ_class_kill(test_grid,test_pft,:,ifm_thin,icut_thin)
                   WRITE(numout,*) 'Clearcut ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_thin,icut_clear)
                   WRITE(numout,*) 'Coppice 1 ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_cop,icut_cop1)
                   WRITE(numout,*) 'Coppice 2 ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_cop,icut_cop2)
                   WRITE(numout,*) 'SRC 1 ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_src,icut_cop1)
                   WRITE(numout,*) 'SRC 2 ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_src,icut_cop2)
                   WRITE(numout,*) 'SRC 3 ',icir,test_pft,&
                        circ_class_kill(test_grid,test_pft,:,ifm_src,icut_cop3)
                ENDIF
             ENDDO

          ENDDO
       ENDDO
    ENDIF
    !-
    

 !! 4. Check mass balance closure
  
    IF (err_act.GT.1) THEN

       ! 3.2 Mass balance closure 
       ! 3.2.1 Calculate final biomass
       pool_end(:,:,:) = zero 
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
          ENDDO
       ENDDO

       ! 3.2.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 3.3 Check mass balance closure
       CALL check_mass_balance("sapiens_forestry_main", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'pft')
      
    ENDIF ! err_act.GT.1


    IF (printlev_loc >= 4) WRITE(numout,*) 'Leaving sapiens_forestry_main'

  END SUBROUTINE sapiens_forestry_main



!! ================================================================================================================================
!! SUBROUTINE    : clearcut_harvest
!!
!>\BRIEF        Marks all the trees at this point and this PFT for clearing
!!
!! DESCRIPTION   : This marks all the trees in this stand for killing.  The
!!                 difference between this and clearcut_noharvest is that
!!                 these trees are moved into the harvest pools in stomate_kill,
!!                 where only the branches should be left onsite.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_kill; the number of trees in each circumference class to kill
!!
!! REFERENCE(S)   : See above, module description.
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{clearcutflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE clearcut_harvest(circ_class_n_temp, circ_class_kill_temp, frac)

    IMPLICIT NONE

 !! 0  Variable and parameter declaration

    !! 0.1 Input variables
    REAL(r_std), DIMENSION(:), INTENT(in)            :: circ_class_n_temp          !! Number of trees in each circumference
                                                                                   !! class
    REAL(r_std), OPTIONAL, INTENT(in)                :: frac                       !! PROGRESSIVE_HARVEST: area fraction of the
                                                                                   !! slot to clearcut (-). Absent = whole slot,
                                                                                   !! i.e. the historical behaviour (bit-neutral).
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(out)           :: circ_class_kill_temp       !! Number of trees within a circ that needs
                                                                                   !! to be killed @tex $(ind m^{-2})$ @endtex

    !! 0.3 Modified variables

    !! 0.4 Local variables


!_ ================================================================================================================================


    ! This routine has gotten much simpler.  All we want to do here is to
    ! indicate how many trees we want to kill.  The actual killing of
    ! the trees will be done in lpj_gap.  As this is a clearcut, we
    ! want to kill all of the trees.
    !
    ! Guillaume M. -- PROGRESSIVE_HARVEST: with ::frac present only that fraction of the
    ! SLOT AREA is clearcut. circ_class_n being a density (ind m-2) over a slot taken as
    ! homogeneous, removing a fraction f of the area removes exactly a fraction f of the
    ! individuals, so scaling the kill IS the area-based clearcut. The matching area
    ! bookkeeping is done downstream in age_class_distr.
    IF (PRESENT(frac)) THEN
       circ_class_kill_temp(:)=frac*circ_class_n_temp(:)
    ELSE
       circ_class_kill_temp(:)=circ_class_n_temp(:)
    ENDIF


  END SUBROUTINE clearcut_harvest

!! ================================================================================================================================
!! SUBROUTINE    : thinning
!!
!>\BRIEF        Thinning of the stand in gridpoint i for PFT j according to the
!! management type (self-thinning, anthropogenic thinning, smoothed anthropogenic
!! thinning or coppicing). Selects which trees are thinned in order to respect 
!! self-thinning rules or rdi_objective. Only trees with circumference lower than 
!! circ_lim can be thinned.
!!
!! DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_ij0; circumference of all individual trees 
!! before thinning for 1 ha at gridpoint i and for PFT j (m), ::vol_thinned 
!! Volume of thinned trees including waste wood \f$(m^3 ha^{-1})\f$ all other readjusted 
!! stand-level variables: ::rdi, ::ind, ::av_height, ...
!!
!! REFERENCE(S)   : See above, module description.
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{thinningflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE thinning(tmin,tmax,thstrat,diameters_temp, &
       circ_temp,circ_class_kill_temp,circ_class_n_temp,&
       rdi_limit,ivm,ipts, ifm)
    
 !! 0  Variable and parameter declaration
    
    !! 0.1 Input variables
    REAL(r_std), INTENT(in)                      :: tmin
    REAL(r_std), INTENT(in)                      :: tmax
    REAL(r_std), INTENT(in)                      :: thstrat
    REAL(r_std), INTENT(in)                      :: rdi_limit
    REAL(r_std), DIMENSION(:), INTENT(in)        :: diameters_temp
    REAL(r_std), DIMENSION(:), INTENT(in)        :: circ_temp
    REAL(r_std), DIMENSION(:), INTENT(in)        :: circ_class_n_temp
    INTEGER(i_std),INTENT(IN)                    :: ivm
    INTEGER(i_std),INTENT(IN)                    :: ifm    
    INTEGER(i_std),INTENT(IN)                    :: ipts


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(out)       :: circ_class_kill_temp

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(ncirc)                :: probability
    REAL(r_std), DIMENSION(ncirc)                :: circ_class_kill_pot
    REAL(r_std)                                  :: circmax
    REAL(r_std)                                  :: circmin
    REAL(r_std)                                  :: current_rdi
    REAL(r_std)                                  :: potential_rdi
    REAL(r_std)                                  :: dead_trees
    REAL(r_std)                                  :: previous_step
    INTEGER    :: icir,endcir,startcir,inccir,nsteps_tried


   INTEGER :: istep1 

!_ ================================================================================================================================    

   ! We are going to assign a thinning probability to every 
   ! circumference class, based on Eq. 12 in Bellassen 2010.  
   ! He does it there for individual trees. The equation 
   ! differs based on the thinning strategy used (from above or below).
   ! Note that although the biomasses are sorted this does not guarantee
   ! that the circumferences are sorted as well because circumference 
   ! largely depends on the heartwood mass but biomass also contains
   ! more volatile pools such as leaves and roots and these pools 
   ! are not always in equilibrium (they may have too little leaves/roots for 
   ! their sapwood). If the model experiences some limitations (for example
   ! C-limitation when doing allocation) there are some rules that first
   ! allocate to a certain circ_class. If this happened the total biomass
   ! of a tree in, e.g. the second circ_class may be higher than total biomass
   ! of a tree in the third circ_class but the circumference of the tree in the
   ! third class will still exceed that of the second class. Long story to
   ! justify why MINVAL and MAXVAL are used instead of circ_temp(1) and 
   ! circ_temp(ncirc)
   circ_class_kill_temp(:)=zero
   circmin=MINVAL(circ_temp(:))
   circmax=MAXVAL(circ_temp(:))
   probability(:)=zero
   
   DO icir=1,ncirc
      
      IF(thstrat .GE. zero)THEN

         ! thinning from below         
         IF(ncirc == 1)THEN

            ! In this case, circmin=circmax and the denominator will be
            ! zero.  However, if we have only one circ class it doesn't really
            ! matter what the probability is set to, because all the biomass
            ! has to come from this circ class.  So we set it equal to the 
            ! minimum probability.
            probability(icir)=tmin

         ELSEIF ((ncirc .GT. 1) .AND. (circmin .EQ. circmax) )THEN

            probability(icir)=tmin
            
         ELSE
            ! Debug
            IF(printlev_loc>=4 .AND. ivm == test_pft .AND. ipts == test_grid)THEN
                    WRITE(numout,*) 'circ_temp(icir): ', circ_temp(icir), &
                    ', icir: ', icir, &
                    ', tmin: ', tmin, &
                    ', tmax: ', tmax, &
                    ', circmin: ', circmin, &
                    ', circmax: ', circmax, &
                    ', thstrat: ', thstrat
                    CALL flush(numout)
            ENDIF

            IF ((circmax-circmin) .LE. min_stomate) THEN
               
               ! There is more than 1 diameter class but the difference between
               ! them is less than min_stomate.
               probability(icir)=tmin
               
            ELSE

               ! There is more than 1 diameter class nd the diameters differ
               probability(icir)=tmin+(tmax-tmin)*((circmax-circ_temp(icir))/&
                    (circmax-circmin))**thstrat

            ENDIF
            
         ENDIF ! ncirc ==1
         
      ELSE !thstrat
         
         ! thinning from above
         IF(ncirc == 1)THEN
          
            ! Same issue as above
            probability(icir)=tmin

         ELSEIF ((ncirc .GT. 1) .AND. (circmin .EQ. circmax) )THEN

            probability(icir)=tmin

         ELSE

            IF ((circmax-circmin) .LE. min_stomate) THEN
               
               ! There is more than 1 diameter class but the difference between
               ! them is less than min_stomate.
               probability(icir)=tmin
               
            ELSE

               ! There is more than 1 diameter class nd the diameters differ
               probability(icir)=tmin+(tmax-tmin)*((circ_temp(icir)-circmin)/&
                    (circmax-circmin))**ABS(thstrat)
               
            ENDIF
            
         ENDIF !ncirc == 1

      ENDIF ! thstrat

      ! Debug
      IF(printlev_loc>=4 .AND. ivm == test_pft .AND. ipts == test_grid)THEN
         WRITE(numout,*) 'probability',probability(icir),icir
         WRITE(numout,*) 'circ_temp',circ_temp(icir),icir
         WRITE(numout,*) 'tmin,tmax',tmin,tmax
         WRITE(numout,*) 'circmin,circmax',circmin,circmax
         WRITE(numout,*) 'thstrat',thstrat
         CALL flush(numout)
      ENDIF
      
      !-
      
   ENDDO

   ! now we need to kill enough trees to get down to the proper density.  If we
   ! are thinning from above, we start from the highest class, and from
   ! below, we start from the lowest.  One thing to keep in mind here is
   ! that it is possible that we have not removed enough trees after looping
   ! through all the circ classes one time, so we need to keep looping through 
   ! until our density is low enough.
   IF(thstrat .GE. zero)THEN
      startcir=1
      endcir=ncirc
      inccir=1
   ELSE
      startcir=ncirc
      endcir=1
      inccir=-1
   ENDIF

   istep1=0
   thin_out:DO
      istep1=istep1+1
      
      IF(istep1 .GT. 1000)THEN

         ! Something is wrong.  Why did it take us so many steps?
         WRITE(numout,*) 'WARNING: Taking too many steps in thinning'
         WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
         WRITE(numout,*) 'current_rdi, rdi_limit: ',current_rdi,rdi_limit
         IF(istep1 .GT. 1001)THEN
            WRITE(numout,*) 'Exiting the loop'
            CALL flush(numout)
            EXIT thin_out
         ENDIF

      ENDIF
      
      DO icir=startcir,endcir,inccir

         ! check the density to see if we have gotten enough trees.  
         ! Notice we only want to include the trees that will still 
         ! be living at the end of this step
         current_rdi=calculate_rdi(diameters_temp(:),circ_class_n_temp(:)-&
                 circ_class_kill_temp(:),alpha_self_thinning(ipts,ivm),ivm,ifm)

         ! Debug
         IF(printlev_loc>=4 .AND. ivm == test_pft .AND. &
              ipts == test_grid)THEN
            WRITE(numout,*) 'Debugging thinning',ivm
            WRITE(numout,*) 'icir,startcir,endcir,inccir',&
                 icir,startcir,endcir,inccir
            WRITE(numout,*) 'rdi_limit',rdi_limit
            WRITE(numout,*) 'current_rdi', current_rdi
            CALL flush(numout)
         ENDIF
         !-

         IF( current_rdi .LE. rdi_limit)THEN      
            EXIT thin_out
         ENDIF

         ! if not, we need to kill some or all of the trees in this class
         ! this is the maximum number of trees we can kill in this circ class
         dead_trees=probability(icir)*&
              (circ_class_n_temp(icir)-circ_class_kill_temp(icir))

         ! if we kill all these trees, can we get the RDI as low as we want?
         circ_class_kill_pot(:)=circ_class_kill_temp(:)
         circ_class_kill_pot(icir)=circ_class_kill_temp(icir)+dead_trees
         potential_rdi=calculate_rdi(diameters_temp(:),&
              circ_class_n_temp(:)-circ_class_kill_pot(:),&
              alpha_self_thinning(ipts,ivm),ivm, ifm)

         ! Debug
         IF(printlev_loc>=4 .AND. ivm == test_pft .AND. ipts == test_grid)THEN
            WRITE(numout,*) 'potential_rdi,rdi_limit',potential_rdi,rdi_limit
            WRITE(numout,*) 'diameters_tmep',diameters_temp(icir)
            WRITE(numout,*) 'circ_class_n_temp',circ_class_n_temp(icir)
            WRITE(numout,*) 'circ_class_kill_pot',circ_class_kill_pot(icir)
            WRITE(numout,*) 'circ_class_kill_temp',circ_class_kill_temp(icir)
            WRITE(numout,*) 'dead_trees',dead_trees
            WRITE(numout,*) 'probability',probability(icir)
            CALL flush(numout)
         ENDIF
         !-

         IF(circ_class_kill_pot(icir) .LT. -min_stomate)THEN
            WRITE(numout,*) 'ERROR: Negative value of circ_class_kill_pot.'
            WRITE(numout,*) 'ipts, ivm, icir: ',ipts,ivm,icir
            WRITE(numout,*) 'pot ',circ_class_kill_pot(icir)
            WRITE(numout,*) 'kill ',circ_class_kill_temp(icir)
            WRITE(numout,*) 'n ',circ_class_n_temp(icir)
            WRITE(numout,*) 'dead_trees ',dead_trees
            WRITE(numout,*) 'probability',probability(icir)
            IF (err_act.GT.1) CALL ipslerr_p (3,'forestry', &
                 'Bad value of circ_class_kill_pot.',&
                 'Look in the output file for ERROR.',&
                 '')
         ENDIF

         IF(potential_rdi .LE. rdi_limit)THEN

            ! yes we can.  How do we know how many trees to kill, though?  
            ! The problem is that the RDI depends on the diameters of the 
            ! remaining trees in a non-linear fashion. Let us do a simple 
            ! bisection search. If we are here, we know that 
            ! circ_class_kill_pot(icir) is too high. So cut it in half and 
            ! try again.
            previous_step=dead_trees/deux
            circ_class_kill_pot(icir)=circ_class_kill_temp(icir)+previous_step
            nsteps_tried=0
            inner_loop:DO

               potential_rdi=calculate_rdi(diameters_temp(:),&
                    circ_class_n_temp(:)-circ_class_kill_pot(:),&
                    alpha_self_thinning(ipts,ivm),ivm, ifm)
               nsteps_tried=nsteps_tried+1
                
               !+++CHECK+++
               ! this value should be externalised
               ! Our RDI is within our tolerance of the target, so we're good.
               IF( ABS(potential_rdi - rdi_limit) .LE. 0.001_r_std)THEN
                  EXIT inner_loop
               ENDIF
               !+++++++++++

               ! Is our potential RDI too high or too low?
               previous_step=previous_step/deux
               IF(potential_rdi .GT. rdi_limit)THEN
                  
                  ! We need to kill more trees, so we add half of our 
                  ! previous step size back.
                  circ_class_kill_pot(icir)=circ_class_kill_pot(icir)+&
                       previous_step
                                
               ELSE
                  
                  ! We need to kill fewer trees, so we take half of our 
                  ! previous step size away.
                  circ_class_kill_pot(icir)=circ_class_kill_pot(icir)-&
                       previous_step

               ENDIF

               ! Debug
               IF(ivm == test_pft .AND. ipts == test_grid .AND. &
                    printlev_loc>=4)THEN
                  WRITE(numout,*) 'nsteps_tried ',nsteps_tried
                  WRITE(numout,*) 'pot ',circ_class_kill_pot(:)
                  WRITE(numout,*) 'kill ',circ_class_kill_temp(:)
                  WRITE(numout,*) 'n ',circ_class_n_temp(:)
                  CALL flush(numout)
               ENDIF
               !-

               !+++CHECK+++
               ! this value should be externalised
               IF(nsteps_tried .GT. 100)THEN
                  WRITE(numout,*) 'ERROR:Trying too many steps in '&
                       //'the bisection search in forestry thinning.'
                  WRITE(numout,*) 'ivm: ',ivm
                  WRITE(numout,*) 'pot ',circ_class_kill_pot(:)
                  WRITE(numout,*) 'kill ',circ_class_kill_temp(:)
                  WRITE(numout,*) 'n ',circ_class_n_temp(:)
                  IF(err_act.GT.1) CALL ipslerr_p (3,'ERROR: forestry', &
                       'Too many steps in the bisection.',&
                       'Look in the output file for ERROR.','')
               ENDIF
               !+++++++++++
               
            ENDDO inner_loop

            ! Schedule the trees for killing.
            circ_class_kill_temp(icir)=circ_class_kill_pot(icir)
            EXIT thin_out

         ELSE

            ! No we can't, so we have to kill all possible trees and 
            ! move onto the next class
            circ_class_kill_temp(icir)=circ_class_kill_temp(icir)+dead_trees
            
         ENDIF

        
      ENDDO
          
      ! check to see if we have any trees left.  We should!  
      ! If not, there is a problem, and we need to catch it 
      ! to avoid an infinite loop.
      IF( (SUM(circ_class_n_temp(:)) - SUM(circ_class_kill_temp(:))) .LE. &
           min_stomate)THEN
         WRITE(numout,*) 'ERROR: We are thinning, but we have no trees left!'
         WRITE(numout,*) 'ivm ',ivm
         WRITE(numout,*) 'kill ',circ_class_kill_temp(:)
         WRITE(numout,*) 'n ',circ_class_n_temp(:)
         IF(err_act.GT.1) CALL ipslerr_p (3,'ERROR: forestry', &
              'Thinning, but no trees left.',&
              'Look in the output file for ERROR.','')
       ENDIF
       
    ENDDO thin_out
   
    ! Debug
    IF(printlev_loc>=4 .AND. ivm == test_pft .AND. ipts == test_grid)THEN
       WRITE(numout,*) 'ivm ',ivm
       WRITE(numout,*) 'circ_kill ',circ_class_kill_temp(:)
       WRITE(numout,*) 'circ_n ',circ_class_n_temp(:)
    END IF
    !-

  END SUBROUTINE thinning
  

!! ================================================================================================================================
!! SUBROUTINE  sapiens_forestry_species_change
!!
!>\BRIEF       When a PFT is harvest a new PFT with another species and management strategy 
!!             is replanted.
!!
!! DESCRIPTION   :   When a PFT is harvest a new PFT with another species and management strategy 
!!                   is replanted. The routine creates new values for veget_max_new which are then
!!                   used in sapiens_lcchange to really change the PFTs.
!!
!!                   The routine sapiens_forestry_species_change is only used when ::ok_change_species is .TRUE. and 
!!                   makes use of the information contained in the species_change map and the 
!!                   desired_fm_map. As an alternative to these maps single values can be specified
!!                   in the run.def for species_change_force and fm_change_force. Although in
!!                   species change and management change are more or less independently called in
!!                   in the code the objective is to call them at the same time. This independence
!!                   was used a bit for debugging but hasn't been tested.
!!                   
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::forest_managed, ::veget_max_new
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================


  SUBROUTINE sapiens_forestry_species_change(npts, lpft_replant, veget_max_new_species, forest_managed, &
       species_change_map, veget_max, fm_change_map)

    IMPLICIT NONE

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                       :: npts                   !! Domain size (unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)       :: species_change_map     !! A map which gives the PFT number that each
                                                                               !! PFT will be replanted as in case of a clearcut.
                                                                               !! (1-nvm,unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)       :: fm_change_map          !! A map which gives the desired FM strategy when
                                                                               !! the PFT will be replanted after a clearcut.
                                                                               !! (1-nvm,unitless)

    !! 0.2 Output variables


    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: veget_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                               !! -> infinity) on ground. Includes
                                                                               !! the nobio fraction, so may sum to
                                                                               !! less than unity.
    REAL(r_std), DIMENSION(:,:),INTENT(out)          :: veget_max_new_species  !! New "maximal" coverage fraction of a PFT 
                                                                               !! (LAI -> infinity)  Includes
                                                                               !! the nobio fraction, so may sum to
                                                                               !! less than unity. (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(inout)           :: lpft_replant           !! Set to true if a PFT has been clearcut
                                                                               !! and needs to be replaced by another species
    INTEGER(i_std), DIMENSION (:,:), INTENT(inout)   :: forest_managed         !! forest management flag (is the forest 
                                                                               !! being managed?)

    !! 0.4 Local variables
    INTEGER                                          :: ipts,ivm,count         !! Indices
    INTEGER                                          :: ivm_young_new_species  !! Indices
    INTEGER                                          :: ivm_young_old_species  !! Indices
    REAL(r_std)                                      :: diff                   !! Difference in surface area before and after 
                                                                               !! a species change (should be zero!)
    REAL(r_std)                                      :: temp_veget_max 
!================================================================================================================================
    
    ! Initialize variable(s)
    veget_max_new_species(:,:)=zero
    count=zero

    ! Check whether replanting is required and if so update
    ! the management strategy where needed.
    DO ipts = 1,npts

       ! Search the PFT that needs to be replanted and replant it
       DO ivm=1,nvm

          IF (is_tree(ivm)) THEN

             ! Was the PFT harvested this year? If so, there is an
             ! opportunity to replant
             IF(lpft_replant(ipts,ivm))THEN

                ! Keep track of the number of species changes. In
                ! case none occur veget_max_new should equal veget_max
                ! else the mass balance check in LCC will think
                ! there was a change in PFTs and will generate errors.
                count = count + 1

                ! Plant the new PFT in the youngest age class. Find
                ! out which PFT is the youngest member of the group 
                ! that we want to regrow.
                ivm_young_new_species = zero
                DO

                   ! Repeat this loop until the matching age class is found
                   ! or an error occurs
                   ivm_young_new_species=ivm_young_new_species+1

                   IF(ivm_young_new_species .GT. nvm)THEN
                      WRITE(numout,*) "ERROR: Can not find young PFT for this age class."
                      WRITE(numout,*) "ipts, ivm: ",ipts, ivm
                      WRITE(numout,*) "species_change_map(ipts,:): ",species_change_map(ipts,:)
                      CALL ipslerr_p (3,'sapiens_forestry_species_change', &
                           'Could not find the youngest PFT in this age class.',&
                           'Something may be wrong with the species_change_map.',&
                           'Check the out_orchidee file for details.')
                   ENDIF

                   ! It is confirmed that we found ivm_young_new_species for the pft
                   ! under consideration. This test requires that the age classes are
                   ! defined such that the youngest age class comes first and that the
                   ! age classes of the same species are grouped (see sapiens_lcchange.f90
                   ! where the age classes are defined)
                   IF(agec_group(ivm_young_new_species) == species_change_map(ipts,ivm))THEN
                      EXIT
                   ENDIF

                ENDDO ! Searching youngest age class for this species

                ! The lpft_replant flag is set in sapiens_forestry or stomate_kill.
                ! Soon after those routines the killing takes place and when
                ! age classes are used veget_max is updated. This means that
                ! the veget_max of the oldest age class that was just harvested
                ! has been moved to the youngest age class of the species that
                ! was harvested. Possibly there was already something in this youngest
                ! age class from previous time steps. The veget_max of the youngest
                ! ages class may thus be larger than the veget_max of the age class 
                ! that was harvested at this time step.
                ivm_young_old_species = start_index(agec_group(ivm))

                ! Is there really a species change?
                IF (ivm .EQ. ivm_young_new_species) THEN

                   ! lpft_replant is true but there is no real species change
                   ! Changes from PFTs earlier in the list could already have
                   ! resulted in an increase of veget_max_new_species
                   !veget_max_new_species(ipts,ivm) = veget_max_new_species(ipts,ivm) + &
                   !     veget_max(ipts,ivm)
                   veget_max_new_species(ipts,ivm_young_new_species) = &
                        veget_max_new_species(ipts,ivm_young_new_species) + &
                        veget_max(ipts,ivm)

                   !---TEMP---
                   IF(printlev_loc>=4) THEN
                      WRITE(numout,*) 'Species change - PFT is the same', ivm
                      WRITE(numout,*) 'Management type, ', &
                           forest_managed(ipts,ivm)
                   ENDIF
                   !----------

                ELSE

                   ! There is a real species change. Add the ground area of the 
                   ! old PFT to the ground area of the new PFT.
                   ! +++CHECK+++
                   ! It seems that in relatively few cases the veget_max of the
                   ! killed PFT has not been moved to the youngest age class
                   ! of that PFT - Need to find when this happens and fix the 
                   ! problem where it occurs (may be it is related to FM3?)
                   ! but here we added a quick patch. Note that
                   ! this code is also used in the case where the youngest age
                   ! class of the PFT was removed. The code moves the whole PFT 
                   ! at once the old PFT can be set to zero.
                   IF (veget_max(ipts,ivm) .GT. min_stomate) THEN

                      ! Replant with the new species
                      veget_max_new_species(ipts,ivm_young_new_species) = &
                           veget_max_new_species(ipts,ivm_young_new_species) + &
                           veget_max(ipts,ivm)
                      veget_max_new_species(ipts,ivm) = zero
                      
                      ! Patch the old species - in land cover change 
                      ! (adjust_delta_veget) we check whether the changes make 
                      ! sense. Those checks will fail unless we make it look
                      ! as if the the youngest PFT was replanted by another
                      ! species
                      IF (ivm .EQ. ivm_young_old_species) THEN
                        WRITE(numout,*) &
                        '1 st class of age cant be replaced by itself' , ivm, ivm_young_old_species
                      ELSE
                        !veget_max(ipts,ivm_young_old_species) = &
                        !    veget_max(ipts,ivm_young_old_species) + veget_max(ipts,ivm)
                        !veget_max(ipts,ivm) = zero
                      ENDIF
                      !---TEMP---
                      IF(printlev_loc>=4) THEN
                         WRITE(numout,*) 'Species change - patch', ivm
                         WRITE(numout,*) 'Future PFT',veget_max(ipts,ivm_young_old_species)
                         WRITE(numout,*) 'current PFT',veget_max(ipts,ivm)
                         WRITE(numout,*) 'Management type, ', forest_managed(ipts,ivm)
                      ENDIF
                      !----------

                   ! Pure case
                   ELSEIF (veget_max(ipts,ivm_young_old_species) .GT. min_stomate) THEN 

                      ! The PFT with the replant label is for example 47. Because the
                      ! PFT has been correctly killed it was replanted in PFT 45 but 
                      ! for PFT 45 the replant flag is FALSE so the surface area will
                      ! be copied from veget_max(ipts,45) to 
                      ! veget_max_new_species(ipts,45). We will first undo the copy
                      ! from 47 to 45. Say that,for example, PFT 45 (young_old) had to 
                      ! become PFT 13 (young_new). We will then replant the correct
                      ! surface are with PFT 13.
                      veget_max_new_species(ipts,ivm_young_old_species) = &
                           veget_max_new_species(ipts,ivm_young_old_species) - &
                           veget_max(ipts,ivm_young_old_species)
                      veget_max_new_species(ipts,ivm_young_new_species) = &
                           veget_max_new_species(ipts,ivm_young_new_species) + &
                           veget_max(ipts,ivm_young_old_species)

                      !---TEMP---
                      IF(printlev_loc>=4) THEN
                         WRITE(numout,*) 'Species change - pure case', ivm
                         WRITE(numout,*) 'management type, ', &
                              forest_managed(ipts,ivm_young_new_species)
                      ENDIF
                      !----------

                   ! Unexpected case
                   ELSE

                      WRITE(numout,*) 'WARNING: species change - unexpected case'
                      IF(err_act.GT.1)THEN
                         WRITE(numout,*) "ipts, ivm: ",ipts, ivm
                         CALL ipslerr_p (3,'sapiens_forestry', &
                           'Species change','unexpected case in IF-loop','')
                      ENDIF

                   ENDIF
                   ! +++++++++++

                ENDIF ! Anthropogenic species change?

                ! Use the new management strategy. The way the desired management
                ! maps were made follows the species logic. A conifer species will
                ! never have FM 3 but a deciduous species can. If a conifer is
                ! changed to a deciduous we need to use the FM of the new deciduous
                ! species in this example FM3. This is reflected in the code by
                ! using the index fom ivm_young_new_species in the left and the right 
                ! hand side. Note that we don't care whether the new PFT has already
                ! some biomass in it or not. Irrespective of its condition we change
                ! FM to the desired strategy.
                forest_managed(ipts,ivm_young_new_species) = &
                     fm_change_map(ipts,ivm_young_new_species)

                !---TEMP---
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'Anthropogenic species change ?'
                   WRITE(numout,*) 'ipts, ivm, ',ipts, ivm  
                   WRITE(numout,*) 'ivm_young_old_species',ivm_young_old_species
                   WRITE(numout,*) 'Species we want to change to, ',species_change_map(ipts,ivm)
                   WRITE(numout,*) 'youngest age class for that species, ', &
                        ivm_young_new_species                  
                   WRITE(numout,*) 'forest_managed, ',forest_managed(ipts,ivm_young_new_species)
                ENDIF
                !----------

             ELSE

                ! No replanting !
                ! veget_max_new will be the same as the old_veget_max.
                ! We add it here just in case another PFT was added to this one
                ! above.
                veget_max_new_species(ipts,ivm) = &
                     veget_max_new_species(ipts,ivm)+veget_max(ipts,ivm)

                ! Given that forest_managed is an inout variable nothing should
                ! be done because we just keep the initial management strategy

             ENDIF

          ELSE

             ! Bare soil, grassland and cropland
             ! veget_max_new will be the same as the old_veget_max.
             ! We add it here just in case another PFT was added to this one
             ! above.
             veget_max_new_species(ipts,ivm) = &
                  veget_max_new_species(ipts,ivm) + veget_max(ipts,ivm)

          ENDIF

       ENDDO !nvm

       ! Debug
       IF(printlev_loc>=4)THEN
          WRITE(numout,*) 'End of routine - final values'
          WRITE(numout,*) 'pixel,PFT,  veget_max, veget_max_new_species'
          DO ivm=1,nvm
             WRITE(numout,*) ipts, ivm,veget_max(ipts,ivm),&
                 veget_max_new_species(ipts,ivm)
          ENDDO

          ! Write warning
          diff = SUM(veget_max_new_species(ipts,:)-veget_max(ipts,:)) 
          IF( ABS(diff) .GT. min_stomate )THEN
             WRITE(numout,*) 'WARNING: species - change surface area is not preserved'
             CALL ipslerr_p (2,'sapiens_forestry', &
                  &          'sapiens_forestry_species_change','surface area is not preserved','')
          ELSE
              WRITE(numout,*) 'Losses and gains cancel each other out'
          ENDIF
          diff = veget_max_new_species(ipts,1)-veget_max(ipts,1)
          IF( ABS(diff) .GT. min_stomate )THEN
             WRITE(numout,*) 'WARNING: species - area of PFT1 is not preserved'
             CALL ipslerr_p (2,'sapiens_forestry', &
                  &          'sapiens_forestry_species_change','PFT1 is not preserved','')
          ELSE
              WRITE(numout,*) 'No leakage to PFT1'
          ENDIF
       ENDIF
       !-

    ENDDO ! npts

    ! Reset this variable
    lpft_replant(:,:)=.FALSE.
    
    ! In case no species changes occured veget_max_new 
    ! should equal veget_max else the mass balance 
    ! check in LCC will think there was a change in 
    ! PFTs and will generate errors.
    IF (count .LT. min_stomate) THEN
       veget_max_new_species(:,:) = veget_max(:,:)
    ENDIF

  END SUBROUTINE sapiens_forestry_species_change


!! ================================================================================================================================
!! SUBROUTINE  : sapiens_forestry_flag_species_change
!!
!>\BRIEF: Check whether there is an opportunity to replant a different species
!!
!! DESCRIPTION : PFTs that died or were harvested should all be empty. This 
!!               is the information that is used to decide whether there is
!!               an opportunity to replant.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
SUBROUTINE sapiens_forestry_flag_species_change(npts, veget_max, circ_class_biomass, circ_class_n, &
     lpft_replant, forest_managed)

 !! 0. Variable and parameter description

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                    :: npts                    !! Number of pixels
    REAL(r_std), DIMENSION(:,:), INTENT(in)       :: veget_max               !! Veget_max at the moment this routine is checked
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in) :: circ_class_biomass      !! ciconference level biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: circ_class_n            !! Density of individuals    
    INTEGER(i_std), DIMENSION (:,:), INTENT(in)   :: forest_managed          !! forest management flag (is the forest 
                                                                             !! being managed?)

    !! 0.2 Output variables

    !! 0.3 Modified variables
    LOGICAL, DIMENSION(:,:), INTENT(inout)        :: lpft_replant            !! Set to true if a PFT has been clearcut
                                                                             !! and needs to be replaced by another species

    !! 0.4 Local variables
    INTEGER(i_std)                                :: ipts,ivm                !! Indices
    INTEGER(i_std)                                :: error                   !! Counter for the numbers of errors
    REAL(r_std)                                   :: biomass
    REAL(r_std)                                   :: ind

!===============================================================================================================================
 
    IF(printlev_loc.GT.4) WRITE(numout,*) 'Entering: change sapiens_forestry_flag_species_change'

    !! 1. Initilaze variables
    error = zero
    !! 2. Check whether there is an opportunity to change the
    !     species. If so flag the opportunity so we can deal 
    !     with the details later.
    DO ivm = 2,nvm

       ! Only forest PFT can be replanted for the moment
       IF(is_tree(ivm))THEN

          DO ipts = 1,npts

            biomass=SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2))
            ind=SUM(circ_class_n(ipts,ivm,:))

             ! If the flag is already true, don't change it
             IF(.NOT.lpft_replant(ipts,ivm))THEN

                ! There is no veget_max. Check whether these PFTs are indeed
                ! empty. If not we have a serious problem
                IF(veget_max(ipts,ivm).LT.min_stomate)THEN

                   IF(biomass.GT.min_stomate .OR. ind.GT.min_stomate)THEN

                      WRITE(numout,*) 'ERROR - no veget_max but biomass and/or n',ipts,ivm
                      WRITE(numout,*) 'biomass, ind, fm, ',biomass,&
                           ind, forest_managed(ipts,ivm)
                      error = error+1

                   ELSEIF(biomass.LT.min_stomate .AND. ind .LT.min_stomate)THEN

                      ! This is the situation we expect. Do nothing
                      
                   ELSE

                      ! May be we overlooked a condition
                      WRITE(numout,*) 'ERROR - overlooked a possible case 1'
                      WRITE(numout,*) 'biomass, n, fm, ',biomass,&
                           ind, forest_managed(ipts,ivm)
                      error = error+1

                   ENDIF
                   
                ELSEIF(veget_max(ipts,ivm).GT.min_stomate)THEN

                   ! There is veget_max. This is where the interesting cases should be
                   IF(biomass.GT.min_stomate .AND. ind.GT.min_stomate)THEN
                   
                      ! This is a simple case: there is still biomass and individuals
                      ! so the stand is alive and there is no opportunity to change
                      ! the species. Note that in sapiens_forestry ::forest_managed
                      ! may already have been changed in case that the desired fm is
                      ! conservation. It is already done in forestry because if 
                      ! fm 2 to 4 is changed to conservation we will not harvest but 
                      ! leave the biomass on site.
                      lpft_replant(ipts,ivm) = .FALSE.
           
                   ELSEIF( biomass.LT.min_stomate .AND. ind.LT.min_stomate)THEN
                      
                      ! The stand was killed or harvest so if the PFT is under
                      ! human management here is an opportunity to replant. 
                      IF(forest_managed(ipts,ivm).EQ.ifm_none)THEN
                      
                         ! If the stand is under a conservation strategy the species will
                         ! not change. Stands under conservation will switch FM without
                         ! a loss of biomass. Hence, there is no chance to replant.
                         ! If the stand is under conservation and dies from natural causes 
                         ! it will be replanted with the same species
                         lpft_replant(ipts,ivm) = .FALSE.

                         ! Debug
                         IF(printlev_loc.GT.4) THEN
                            WRITE(numout,*) 'Replant - no biomass & no ind, conservation',ipts,ivm
                            WRITE(numout,*) 'biomass, ind, fm, ',biomass,ind, forest_managed(ipts,ivm)
                            WRITE(numout,*) 'lpft_replant, ',lpft_replant(ipts,ivm)
                         END IF
                         ! -

                      ELSEIF(forest_managed(ipts,ivm).NE.ifm_none)THEN
                      
                         ! The stand was managed. This is a chance to change the species.
                         ! May be we wont use the chance but that will only be checked at
                         ! the end of the year in sapiens_forestry_species_change (sapiens_forestry.f90) 
                         lpft_replant(ipts,ivm) = .TRUE.

                         ! Debug
                         IF(printlev_loc.GT.4) THEN
                            WRITE(numout,*) 'Replant - no biomass & no ind, management',ipts,ivm
                            WRITE(numout,*) 'biomass, ind, fm, ',biomass,ind, forest_managed(ipts,ivm)
                            WRITE(numout,*) 'lpft_replant, ',lpft_replant(ipts,ivm)

                         END IF
                         !-

                      ELSE

                         WRITE(numout,*) 'ERROR - overlooked a possible option 1', ipts, ivm
                         WRITE(numout,*) 'biomass, ind, fm, ',biomass, ind, forest_managed(ipts,ivm)
                         WRITE(numout,*) 'lpft_replant, ',lpft_replant(ipts,ivm)
                         error = error+1

                      ENDIF

                   ELSE

                      WRITE(numout,*) 'ERROR - overlooked a possible case 2'
                      WRITE(numout,*) 'biomass, ind, fm, ',biomass,ind, forest_managed(ipts,ivm)
                      WRITE(numout,*) 'lpft_replant, ',lpft_replant(ipts,ivm)
                      error = error+1
                   ENDIF

                ELSE

                   WRITE(numout,*) 'ERROR - overlooked a possible case 3'
                   WRITE(numout,*) 'biomass, ind, fm, ',biomass,ind, forest_managed(ipts,ivm)
                   WRITE(numout,*) 'lpft_replant, ',lpft_replant(ipts,ivm)
                   error = error+1
                   
                ENDIF

             ENDIF ! lpft_replant
          
          ENDDO ! ivm

       ENDIF ! is_tree
       
    ENDDO ! ipts

    IF(error.GT.zero)THEN
       CALL ipslerr_p(3,'sapiens_forestry_flag_species_change','ERROR: this subroutine did not function as expected','','')
    ENDIF

  END SUBROUTINE sapiens_forestry_flag_species_change


!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_set_fm
!!
!>\BRIEF        Set forest managment strategy for each pixel and each PFT.
!!              This is done by reading a map or by using parameters from run.def
!!
!! DESCRIPTION   :  
!!     
!!
!! FM = 1 : Unmanaged, no human intervention (ORCHIDEE default)
!!      2 : Thinnings based on the RDI, clearcuts based on tree density,
!!              annual increment, and tree diameter.  Thinnings from above and
!!              from below are determined by the sign of thstrat.
!!      3 : Coppices
!!      4 : Short rotation coppices
!!      5 : Continuous cover forestry

!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::forest_managed
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_set_fm ( npts, lalo, neighbours, resolution, contfrac, forest_managed )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(out)    :: forest_managed  !! Forest management flag: 1 = unmanaged (orchidee 
                                                                           !! standard), 2= self-thinning only, 3= coppices,
                                                                           !! 4= short rotation coppices, 5= continuous
                                                                           !! cover forestry
    !! 0.3 Local variables
    INTEGER(i_std)                                      :: ivm, ipts
    INTEGER(i_std), DIMENSION(npts,nvmap)               :: fm_map_temp     !! A temporary variable to hold the forest

    ! Initialize printlev_loc
    IF (firstcall_sapiens_forestry) THEN
       ! Initialize local printlev if it is not already done
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_forestry=.FALSE.
    END IF

    
    ! Set forest_managed depending on the option ok_read_fm_map: read from file or from run.def
    IF (ok_read_fm_map) THEN
       ! forest_managed will be read from file.
       IF (printlev_loc >= 3) WRITE(numout,*) 'In sapiens_forestry_set_fm, set forest_managed by reading a map'
       
       ! If we are using age classes, we read in the map in the same
       ! way but then we change it a bit to account for age classes.
       CALL sapiens_forestry_read_fm(npts, lalo, neighbours, resolution, contfrac, fm_map_temp)
       
       IF (nagec .GT. 1) THEN
          ! All age classes of the same PFT will have the same management
          forest_managed(:,:) = zero
          DO ivm = 1,nvm
             forest_managed(:,ivm)=fm_map_temp(:,agec_group(ivm))
          ENDDO
       ELSE
          forest_managed(:,:)=fm_map_temp(:,:)
       ENDIF

       !+++CHECK+++
       ! The global forest management reconstructions distinguish unmanaged and
       ! managed forests. For the moment the only management in these maps is a
       ! rotational even aged management with thinnings and fellings. At the site
       ! level a thinning easily removes 30% of the biomass and a felling removes 
       ! 100%. If the same is applied at large pixels (2500 to 40000 km) this 
       ! results is unrealistic biomass removals with consequences for the albedo
       ! transpiration and roughness of the entire pixel. When all the hectare-scale
       ! stands in a pixel are managed, we think the pixel-level management would
       ! look like continuous cover forestry in the sense that the pixel as a whole
       ! is never harvested entirely in a single year. As long as we are not 100%
       ! sure of this solution yet, we keep the original management maps but 
       ! reinterpret ifm_thin (thin and fell or even-aged rotational management) as 
       ! ifm_uneven (continuous cover forestry).
       IF (concept_scale=='global') THEN
          WHERE(forest_managed==ifm_thin)
             forest_managed(:,:) = ifm_uneven
          END WHERE
       END IF
       !+++++++++++

       ! Guillaume M. -- SHORT ROTATION COPPICE: PFTs flagged FM_SRC_PFT switch to ifm_src
       ! wherever the map says MANAGED; unmanaged cells (ifm_none) are left alone, since a
       ! coppice implies an intervention. Eucalyptus is managed by stump resprouting, which
       ! is neither a plantation nor natural regeneration in the model's sense.
       !
       ! Guillaume M. -- /!\ ifm_src was already FULLY implemented (*_src parameters,
       ! icut_cop1/2/3, first/second_coppice_cut) but NEVER assigned, for want of an
       ! assignment rule. /!\ Must stay AFTER the global reinterpretation above, which would
       ! otherwise overwrite it by turning every ifm_thin back into ifm_uneven.
       IF (ALLOCATED(fm_src_pft)) THEN
          DO ivm = 1, nvm
             IF (.NOT. fm_src_pft(ivm)) CYCLE
             WHERE (forest_managed(:,ivm) /= ifm_none)
                forest_managed(:,ivm) = ifm_src
             END WHERE
             IF (printlev_loc >= 2) WRITE(numout,*) &
                  'sapiens_forestry_set_fm : PFT ', ivm, ' passe en taillis (ifm_src), ', &
                  COUNT(forest_managed(:,ivm) == ifm_src), ' mailles'
          END DO
       END IF

       ! Guillaume M. -- Same rule as FM_SRC_PFT above, for the LONG rotation coppice
       ! (ifm_cop). Holm oak is managed as a 35-60 yr simple coppice, not as an even-aged
       ! high forest. Unmanaged cells are left alone: a coppice implies an intervention.
       ! Must stay AFTER the global ifm_thin -> ifm_uneven reinterpretation, which would
       ! otherwise overwrite it.
       IF (ALLOCATED(fm_cop_pft)) THEN
          DO ivm = 1, nvm
             IF (.NOT. fm_cop_pft(ivm)) CYCLE
             WHERE (forest_managed(:,ivm) /= ifm_none)
                forest_managed(:,ivm) = ifm_cop
             END WHERE
             IF (printlev_loc >= 2) WRITE(numout,*) &
                  'sapiens_forestry_set_fm : PFT ', ivm, ' passe en taillis simple (ifm_cop), ', &
                  COUNT(forest_managed(:,ivm) == ifm_cop), ' mailles'
          END DO
       END IF
       
    ELSE
       
       ! Initialize using values for FOREST_MANAGED_FORCED from run.def
       IF (printlev_loc >= 3) WRITE(numout,*) 'In sapiens_forestry_set_fm, set forest_managed with parameters from run.def'
       DO ipts = 1, npts
          forest_managed(ipts,:) = forest_managed_forced(:)
       END DO
    END IF
          
  END SUBROUTINE sapiens_forestry_set_fm
  
!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_read_fm
!!
!>\BRIEF        Read in a map that gives the forest management strategy
!!              for each pixel and each PFT.
!!
!! DESCRIPTION   :  
!!     
!!
!! FM = 1 : Unmanaged, no human intervention (ORCHIDEE default)
!!      2 : Thinnings based on the RDI, clearcuts based on tree density,
!!              annual increment, and tree diameter.  Thinnings from above and
!!              from below are determined by the sign of thstrat.
!!      3 : Coppices
!!      4 : Short rotation coppices
!!      5 : Continuous cover forestry
!!
!!       NOTE: This routine was mostly copied from slowproc where the PFTmap is read in.
!!             Grid interpolation is used, but only to look at the nearby pixels to see
!!             see which management strategy is dominant.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::forest_managed
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_read_fm ( npts, lalo, neighbours, resolution, contfrac, forest_managed )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)         :: forest_managed  !! Forest management flag: 1 = unmanaged (orchidee 
                                                                           !! standard), 2= self-thinning only, 3= coppices,
                                                                           !! 4= short rotation coppices, 5= continous
                                                                           !! cover forestry

    !! 0.3 Modified fields


    !! 0.4 Local variables
    LOGICAL                                                :: debug=.FALSE.   !! A flag to print out debugging information.
    INTEGER(i_std)                                         :: fid             !! The ID of the NetCDF file.
    INTEGER(i_std)                                         :: nb_coord        !! The number of coordinates in the NetCDF file    
    INTEGER(i_std)                                         :: nb_gat          !!
    INTEGER(i_std)                                         :: nb_var          !! The number of variables in the NetCDF file
    INTEGER(i_std)                                         :: nb_dim          !! The number of dimensions in the NetCDF file
    INTEGER(i_std)                                         :: iml, jml, lml,ivm !! indices
    INTEGER(i_std)                                         :: ip, inbv, jp    !! indices
    LOGICAL                                                :: l_ex            !! A flag which indicates if a variable
                                                                              !! exists in the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                 :: lat_lu, lon_lu  !! The latitude and longitude read in from
                                                                              !! the NetCDF file    
    INTEGER,DIMENSION(flio_max_var_dims)                   :: l_d_w           !! List of the dimension lengths of the variable
                                                                              !! in the NetCDF file
    INTEGER(i_std)                                         :: ipts            !! index
    INTEGER(i_std)                                         :: alloc_err       !! A flag tripped if we have an error in allocation
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE               :: fmmap_r         !! The map read in from the NetCDF file
    INTEGER(i_std),DIMENSION(:,:,:),ALLOCATABLE            :: fmmap_i         !! The integer form of the map read in
    INTEGER(i_std)                                         :: closest_lat     !! The index of the closest latitude we found.
    INTEGER(i_std)                                         :: closest_lon     !! The index of the closest longitude we found.
    REAL(r_std)                                            :: distance        !! The distance from the current point to the 
                                                                              !! point on the map
    REAL(r_std)                                            :: closest_dist    !! The distance to the closet point we've found.
    INTEGER                                                :: large_int       !! A number which indicates that the grid data
                                                                              !! is not available
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_ful, lon_ful
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: ibvm, nbvmax,ifm, ivma
    LOGICAL ::           ok_interpol ! optionnal return of aggregate_2d
    REAL(r_std)                                            :: sum_fm,sumf
    REAL(r_std),DIMENSION(nfm_types)                       :: fm_sum
    REAL(r_std), DIMENSION(npts,nvmap,nfm_types)           :: fm_frac ! (nfm_types=7)
    LOGICAL                                                :: ltemp ! temporary logical variable
    INTEGER(i_std)                                         :: fm_default    !! Default value used if no information was found in the input map
    !_ ================================================================================================================================

    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Entering sapiens_forestry_read_fm'

    ! This integer has to be large enough that it never shows up on the map, without being
    ! so large that is causes overflows.  Since none of the points on the map should be
    ! larger than the number of FM strategies we have (nfm_types), this is sufficiently big.
    large_int = nvmap*nfm_types+1

    ! Set 1 as default value to use if there is no information in the input file for a specific pixel.
    fm_default=1

    IF (xios_orchidee_ok .AND. xios_interpolation_fm) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "sapiens_forestry_read_fm: Use XIOS to read and interpolate " &
            // TRIM(filename_fm) // " for variable FM_STRAT"

       ! Read the fraction from the input file for each of the forest management methodes. 
       ! fm_frac1 is the fraction of covered by the method 1, etc.
       CALL xios_orchidee_recv_field('fm_frac1',fm_frac(:,:,1))
       CALL xios_orchidee_recv_field('fm_frac2',fm_frac(:,:,2))
       CALL xios_orchidee_recv_field('fm_frac3',fm_frac(:,:,3))
       CALL xios_orchidee_recv_field('fm_frac4',fm_frac(:,:,4))
       CALL xios_orchidee_recv_field('fm_frac5',fm_frac(:,:,5))
       CALL xios_orchidee_recv_field('fm_frac6',fm_frac(:,:,6))
       CALL xios_orchidee_recv_field('fm_frac7',fm_frac(:,:,7))
       
       DO ipts = 1, npts
          ! We need to do this for every PFT from file
          DO ivma=1,nvmap
             ! If it's not a forest, we don't care.
             IF(.NOT. is_tree(start_index(ivma)))THEN
                ! Current pft is not a tree so we don't care, set zero
                forest_managed(ipts,ivma)=0
             ELSE IF (SUM(fm_frac(ipts,ivma,:)) <= min_sechiba) THEN
                ! This is a tree pft but there is no method in the input file
                forest_managed(ipts,ivma)=fm_default
             ELSE
                ! Take the index which has the biggest fraction
                ! The index corresponds to the forest management method red from file
                forest_managed(ipts,ivma)=MAXLOC(fm_frac(ipts,ivma,:),1)
             END IF
          END DO
       END DO

       IF (printlev>=2) WRITE(numout,*) 'Done with interpolating the FM map using XIOS'
       
    ELSE

       ! Read and interpoaltion using IOIPSL and aggregate_p (without XIOS)
       ! Get file name again because sapiens_forestry_xios_initialize might not be called (installation without XIOS)
       filename_fm = 'FMmap.nc'
       CALL getin_p('FM_FILE',filename_fm)

       IF (printlev_loc >= 1) WRITE(numout,*) "sapiens_forestry_read_fm: Use IOIPSL to read and aggregate_p to interpolate " &
            // TRIM(filename_fm) // " for variable FM_STRAT"
       IF (is_root_prc) THEN
          IF (debug) THEN
             WRITE(numout,*) "Entering sapiens_forestry_read_fm. Debug mode."
             WRITE (*,'(/," --> fliodmpf")')
             CALL fliodmpf (TRIM(filename_fm))
             WRITE (*,'(/," --> flioopfd")')
          ENDIF
          CALL flioopfd (TRIM(filename_fm),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
          IF (debug) THEN
             WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
             WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
             WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
          ENDIF
       ENDIF
       CALL bcast(nb_coord)
       CALL bcast(nb_var)
       CALL bcast(nb_gat)
       
       ! This finds the number of longitude points in the file.
       IF (is_root_prc) &
            CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
       CALL bcast(l_d_w)
       iml=l_d_w(1)
       WRITE(numout,*) "FM Map: iml =",iml
       
       ! This finds the number of latitude points in the file.
       IF (is_root_prc) &
            CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
       CALL bcast(l_d_w)
       jml=l_d_w(1)
       WRITE(numout,*) "FM Map: jml =",jml
       
       ! Now find the number of PFTs in the file.  If this is not equal to the number
       ! of PFTs that we actually have, that's a problem and we'll crash.
       IF (is_root_prc) &
            CALL flioinqv (fid,v_n="FM_STRAT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
       CALL bcast(l_d_w)
       lml=l_d_w(3)
       
       IF (lml /= nvmap) THEN
          WRITE(numout,*) 'lml = ',lml
          WRITE(numout,*) 'nvmap = ',nvmap
          WRITE(numout,*) 'Stopping. '
          CALL ipslerr_p (3,'sapiens_forestry', &
               'Problem with forest management strategy map.','lml /= nvmap', &
               '(number of pft must be equal)')
       ENDIF
       
       ! Allocate the map that will be read in
       WRITE(numout,*) 'Reading the Forest Management strategy file'
       
       ALLOCATE(fmmap_r(iml,jml,nvmap), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
         'Problem in allocation of variable fmmap_r','','')

       ALLOCATE(fmmap_i(iml,jml,nvmap), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
            'Problem in allocation of variable fmmap_i','','')
      
       ! This reads in the map that is in the file
       IF (is_root_prc) THEN
          fmmap_r(:,:,:)=large_int*2.0
          CALL fliogetv (fid,"FM_STRAT", fmmap_r, start=(/ 1, 1, 1 /), count=(/ iml, jml, nvmap /))
          ! Right now the values are all real, but they should be integers.
          ! Careful, NINT might not work if the precision of fmmap_r is not single.
          ! In that case, IDNINT should be used.
          DO ip=1,iml
             DO jp=1,jml
                DO ivma=1,nvmap
                   ! There will be some fill values in here.  If we pass
                   ! a large fill value to NINT, it crashes.  So let's 
                   ! test for it
#ifdef CPP_IEEE_ARITHMETIC
                   ltemp=IEEE_IS_NAN(fmmap_r(ip,jp,ivma))    
#else
                   ltemp=isnan(fmmap_r(ip,jp,ivma))    
#endif
                   IF(ltemp)THEN
                      fmmap_i(ip,jp,ivma)=large_int
                   ELSE
                      IF(fmmap_r(ip,jp,ivma) .GE. 0.0 .AND. fmmap_r(ip,jp,ivma) < large_int)THEN
                         fmmap_i(ip,jp,ivma)=NINT(fmmap_r(ip,jp,ivma))
                      ELSE
                         ! This value should be big enough that we don't barl ourselves
                         ! below.
                         fmmap_i(ip,jp,ivma)=large_int
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       
       CALL bcast(fmmap_i)
       
       ! Now I need to get the latitude and longitude
       ! First, get the axes from the map file.
       ALLOCATE(lat_lu(jml), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
         'Problem in allocation of variable lat_lu','','')

       ALLOCATE(lon_lu(iml), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
         'Problem in allocation of variable lon_lu','','')

       IF (is_root_prc) THEN
          CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
       ENDIF
       CALL bcast(lon_lu)
       CALL bcast(lat_lu)
       
       ! Now I can interpolate.  Remember that we are dealing with integer
       ! values and management strategies.  Therefore, we cannot do a strict
       ! interpolation, since that might gives values of the strategy to be
       ! 2.3, or something ridiculous.  We cannot just round the number, either,
       ! since we could have some squares with 1 (no management) and some with
       ! 3 (coppices);  An average of that would be 2 (high stands), which doesn't
       ! make any sense.  So we look to see what the most prevalent type of
       ! management of nearby pixels is, and then just take that.
       ALLOCATE(lat_ful(iml,jml), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
            'Problem in allocation of variable lat_ful','','')

       ALLOCATE(lon_ful(iml,jml), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
            'Problem in allocation of variable lon_ful','','')

       DO ip=1,iml
          lon_ful(ip,:)=lon_lu(ip)
       ENDDO
       DO jp=1,jml
          lat_ful(:,jp)=lat_lu(jp)
       ENDDO
       !
       ! Mask of permitted variables.
       !
       ALLOCATE(mask(iml,jml), STAT=alloc_err)
       IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
            'Problem in allocation of variable mask','','')

       mask(:,:) = 0
       DO ip=1,iml
          DO jp=1,jml
             ! If at least one PFT has a management strategy, we should interpolate
             ! this grid point.
             sum_fm=SUM(fmmap_i(ip,jp,:))
             IF ( sum_fm .GT. min_sechiba .AND. sum_fm .LT. large_int) THEN
                mask(ip,jp) = 1
                IF (debug) THEN
                   WRITE(numout,*) "update : SUM(fmmap(",ip,jp,")) = ",sum_fm
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       
       !
       ! The number of maximum vegetation map points in the GCM grid should
       ! also be computed and not imposed here.
       !
       nbvmax = 200
       
       callsign="Forest Management map"
       
       ok_interpol = .FALSE.
       DO WHILE ( .NOT. ok_interpol )
          WRITE(numout,*) "Projection arrays for ",callsign," : "
          WRITE(numout,*) "nbvmax = ",nbvmax

          ALLOCATE(sub_index(npts, nbvmax,2), STAT=alloc_err)
          IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
               'Problem in allocation of variable sub_index','','')
          sub_index(:,:,:)=0
          
          ALLOCATE(sub_area(npts, nbvmax), STAT=alloc_err)
          IF (alloc_err /= 0) CALL ipslerr_p(3,'sapiens_forestry_read_fm',&
               'Problem in allocation of variable sub_area','','')
          sub_area(:,:)=zero
          
          CALL aggregate_p(npts, lalo, neighbours, resolution, contfrac, &
               iml, jml, lon_ful, lat_ful, mask, callsign, &
               nbvmax, sub_index, sub_area, ok_interpol)
          
          IF ( .NOT. ok_interpol ) THEN
             DEALLOCATE(sub_area)
             DEALLOCATE(sub_index)
          ENDIF
          
          nbvmax = nbvmax * 2
       ENDDO

       ! The FM maps with 28 PFTs for Europe were made on the exact grid that
       ! is used in the FG4 simulations and therefore avoids the need for 
       ! interpolation. Moreover, the background value (used outside Europe and
       ! within Europe where the PFT is not present) was no fm (=1). Hence, ifm
       ! can never by zero (as it should be). Because all difficulties were
       ! tackled when making the FM maps, reading the maps is straightforward.
       ! The FM maps with 15 PFTs for the world were made on a 0.25 grid which
       ! has a much higher resolution than most of the applications (FG1, 
       ! FG2, and FG3), hence, regridding will be necessary. Also to make better
       ! looking maps land without the specific forest PFT got FM strategy zero. 
       ! Regridding comes with at least two risks: (1) the FM maps use the most 
       ! frequent FM strategy within the search area. This could be "no 
       ! management" because the minority of the pixels in the search area 
       ! contains forest. At the same time the PFT maps will also need to be 
       ! regridded but they use an average value. This way it seems possible that 
       ! we have inconsistencies between the PFT maps (forest present) and the FM 
       ! maps (no management strategy available -> fm = zero). (2) We could also
       ! have used 1 as the background values. The maps would then look less 
       ! realistic but for the simulations itself it would make little difference
       ! because the PFT is not present at the pixels were the background values
       ! are set. When regridding this could result in an edge effect of having 
       ! some unmanaged pixels at the edge of a region with  managed forests. The
       ! code below is suitable for reading both types of FM maps (FG4 vs 
       ! FG1, FG2 and FG3).
       DO ipts = 1, npts
          ! For this point, we need to see what dominant FM types are
          ! nearby.
          sumf=zero
          ! We need to do this for every PFT
          DO ivma=1,nvmap
             ! If it's not a forest, we don't care.
             IF(.NOT. is_tree(start_index(ivma)))THEN
                forest_managed(ipts,ivma)=0
             ELSE
                fm_sum(:)=0.0
                DO ibvm=1, nbvmax
                   ! Leave the do loop if all sub areas are treated, sub_area <= 0
                   IF ( sub_area(ipts,ibvm) <= zero ) EXIT
                   ip = sub_index(ipts,ibvm,1)
                   jp = sub_index(ipts,ibvm,2)
                   ifm=fmmap_i(ip,jp,ivma)
                   ! Only calculate the surface area of the forested pixels.
                   ! This approach should ensure that the regridding of the
                   ! FM maps and the PFT maps will use exactly the same grid.
                   IF (ifm.NE.zero) THEN
                      fm_sum(ifm)=fm_sum(ifm)+sub_area(ipts,ibvm)
                   END IF
                ENDDO
                ! Whichever FM type has the most area, we use that one.
                IF (SUM(fm_sum(:)).EQ.zero) THEN
                   ! This is a pixel without forest. Because of the consistency
                   ! between the PFT and the FM maps it is safe to set FM to zero.
                   ! The biggest advantage of doing so is that the variable
                   ! forest_managed should show an intuitive map when being plotted.
                   forest_managed(ipts,ivma) = fm_default
                ELSE
                   forest_managed(ipts,ivma)=MAXLOC(fm_sum(:),1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       
       IF (printlev>=2) WRITE(numout,*) 'Done with interpolating the FM map with IOIPSL and aggregate_p method'
       
       
       DEALLOCATE(fmmap_i)
       DEALLOCATE(fmmap_r)
       DEALLOCATE(lat_lu,lon_lu)
       DEALLOCATE(lat_ful,lon_ful)
       DEALLOCATE(mask)
       IF(ALLOCATED(sub_index)) DEALLOCATE(sub_index)
       IF(ALLOCATED(sub_area)) DEALLOCATE(sub_area)
       
    END IF ! xios_interpolation

    CALL xios_orchidee_send_field("interp_diag_forest_managed",REAL(forest_managed)) 

    IF ( printlev_loc >= 5 ) WRITE(numout,*) 'Leaving sapiens_forestry_read_fm'

 END SUBROUTINE sapiens_forestry_read_fm

!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_read_spinup_clearcut
!!
!>\BRIEF        Read in a map that gives the state of clearcut during spinup
!!              for each pixel and each PFT.
!!
!! DESCRIPTION   :  
!!     
!!
!!
!!       NOTE: This routine was mostly copied from slowproc where the PFTmap is read in.
!!             Grid interpolation is used, but only to look at the nearby pixels to see
!!             see which management strategy is dominant.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::spinup_clearcut
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_read_spinup_clearcut ( npts, lalo, neighbours, resolution, contfrac, spinup_clearcut )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)         :: spinup_clearcut  !! Forest management flag: 1 = unmanaged (orchidee
                                                                           !! standard), 2= self-thinning only, 3= coppices,
                                                                           !! 4= short rotation coppices, 5= continous
                                                                           !! cover forestry
    !! 0.3 Modified fields


    !! 0.4 Local variables
    CHARACTER(LEN=80)                                      :: filename        !! A string to hold the file name
    LOGICAL                                                :: debug=.FALSE.   !! A flag to print out debugging information.
    INTEGER(i_std)                                         :: fid             !! The ID of the NetCDF file.
    INTEGER(i_std)                                         :: nb_coord        !! The number of coordinates in the NetCDF file    
    INTEGER(i_std)                                         :: nb_gat          !!
    INTEGER(i_std)                                         :: nb_var          !! The number of variables in the NetCDF file
    INTEGER(i_std)                                         :: nb_dim          !! The number of dimensions in the NetCDF file
    INTEGER(i_std)                                         :: iml, jml, lml,ivm !! indices
    INTEGER(i_std)                                         :: ip, inbv, jp    !! indices
    LOGICAL                                                :: l_ex            !! A flag which indicates if a variable
                                                                              !! exists in the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                 :: lat_lu, lon_lu  !! The latitude and longitude read in from
                                                                              !! the NetCDF file    
    INTEGER,DIMENSION(flio_max_var_dims)                   :: l_d_w           !! List of the dimension lengths of the variable
                                                                              !! in the NetCDF file
    INTEGER(i_std)                                         :: ipts            !! index
    INTEGER(i_std)                                         :: ALLOC_ERR       !! A flag tripped if we have an error in allocation
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE               :: fmmap_r         !! The map read in from the NetCDF file
    INTEGER(i_std),DIMENSION(:,:,:),ALLOCATABLE            :: fmmap_i         !! The integer form of the map read in
    INTEGER(i_std)                                         :: closest_lat     !! The index of the closest latitude we found.
    INTEGER(i_std)                                         :: closest_lon     !! The index of the closest longitude we found.
    REAL(r_std)                                            :: distance        !! The distance from the current point to the 
                                                                              !! point on the map
    REAL(r_std)                                            :: closest_dist    !! The distance to the closet point we've found.
    INTEGER                                                :: large_int       !! A number which indicates that the grid data
                                                                              !! is not available
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_ful, lon_ful
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: ibvm, nbvmax,ifm, ivma
    LOGICAL ::           ok_interpol ! optionnal return of aggregate_2d
    REAL(r_std)                                            :: sum_fm,sumf
    LOGICAL ::           ltemp ! temporary logical variable
  !_ ================================================================================================================================

    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Entering sapiens_forestry_read_spinup_clearcut'

    ! This integer has to be large enough that it never shows up on the map, without being
    ! so large that is causes overflows.  Since none of the points on the map should be
    ! larger than the number of FM strategies we have (nfm_types), this is sufficiently big.
    large_int = flag_spinup_clearcut + 1

   !
    !Config Key   = SPINUP_CLEARCUT_FILE
    !Config Desc  = Name of file to be read
    !Config If    = OK_STOMATE
    !Config Def   = spinup_clearcut.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'spinup_clearcut.nc'
    CALL getin_p('SPINUP_CLEARCUT_FILE',filename)

    IF (is_root_prc) THEN
       IF (debug) THEN
          WRITE(numout,*) "Entering sapiens_forestry_read_spinup_clearcut. Debug mode."
          WRITE (*,'(/," --> fliodmpf")')
          CALL fliodmpf (TRIM(filename))
          WRITE (*,'(/," --> flioopfd")')
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (debug) THEN
          WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF
    CALL bcast(nb_coord)
    CALL bcast(nb_var)
    CALL bcast(nb_gat)

    ! This finds the number of longitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    iml=l_d_w(1)
    WRITE(numout,*) "spinup_clearcut map: iml =",iml

    ! This finds the number of latitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    jml=l_d_w(1)
    WRITE(numout,*) "spinup_clearcut map: jml =",jml

    ! Now find the number of PFTs in the file.  If this is not equal to the number
    ! of PFTs that we actually have, that's a problem and we'll crash.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="Clearcut",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    lml=l_d_w(3)

    IF (lml /= nvmap) THEN
       WRITE(numout,*) 'lml = ',lml
       WRITE(numout,*) 'nvmap = ',nvmap
       WRITE(numout,*) 'Stopping. '
       CALL ipslerr_p (3,'sapiens_forestry', &
            &          'Problem with spinup clearcut map.','lml /= nvmap', &
            &          '(number of pft must be equal)')
    ENDIF
    !

    ! Allocate the map that will be read in
    WRITE(numout,*) 'Reading the spinup clearcut strategy file'
    !
    ALLOC_ERR=-1
    ALLOCATE(fmmap_r(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of fmmap_r : ",ALLOC_ERR
       
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(fmmap_i(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of fmmap_i : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in fmmpap_i','','')
    ENDIF

    ! This reads in the map that is in the file
    IF (is_root_prc) THEN
       fmmap_r(:,:,:)=zero
       CALL fliogetv (fid,"Clearcut", fmmap_r, start=(/ 1, 1, 1 /), count=(/ iml, jml, nvmap /))
       ! Right now the values are all real, but they should be integers.
       ! Careful, NINT might not work if the precision of fmmap_r is not single.
       ! In that case, IDNINT should be used.
       DO ip=1,iml
          DO jp=1,jml
             DO ivma=1,nvmap
                ! There will be some fill values in here.  If we pass
                ! a large fill value to NINT, it crashes.  So let's 
                ! test for it
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(fmmap_r(ip,jp,ivma))    
#else
                ltemp=isnan(fmmap_r(ip,jp,ivma))    
#endif
                IF(ltemp)THEN
                   fmmap_i(ip,jp,ivma)=zero
                ELSE
                   IF(fmmap_r(ip,jp,ivma) .GE. 0.0 .AND. fmmap_r(ip,jp,ivma) .LE. flag_spinup_clearcut)THEN
                      fmmap_i(ip,jp,ivma)=NINT(fmmap_r(ip,jp,ivma))
                   ELSE
                      ! For the moment there is a strict rule that values of input 
                      ! cannot be bigger than 1.
                      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'value is bigger than 1','','') 
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL bcast(fmmap_i)

    ! Now I need to get the latitude and longitude
    ! First, get the axes from the map file.
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in lat_lu','','') 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in lon_lu','','') 
    ENDIF
    IF (is_root_prc) THEN
       CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)

    ! Now I can interpolate.  Remember that we are dealing with integer
    ! values and management strategies.  Therefore, we cannot do a strict
    ! interpolation, since that might gives values of the strategy to be
    ! 2.3, or something ridiculous.  We cannot just round the number, either,
    ! since we could have some squares with 1 (no management) and some with
    ! 3 (coppices);  An average of that would be 2 (high stands), which doesn't
    ! make any sense.  So we look to see what the most prevalent type of
    ! management of nearby pixels is, and then just take that.
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in lat_ful','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in lon_ful','','')
    ENDIF
    !
    DO ip=1,iml
       lon_ful(ip,:)=lon_lu(ip)
    ENDDO
    DO jp=1,jml
       lat_ful(:,jp)=lat_lu(jp)
    ENDDO
    !
    ! Mask of permitted variables.
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in mask','','')
    ENDIF
    !
    mask(:,:) = 1
    ! DO ip=1,iml
    !    DO jp=1,jml
    !       ! If at least one PFT has a management strategy, we should interpolate
    !       ! this grid point.
    !       sum_fm=SUM(fmmap_i(ip,jp,:))
    !       IF ( sum_fm .GT. min_sechiba .AND. sum_fm .LT. large_int) THEN
    !          mask(ip,jp) = 1
    !          IF (debug) THEN
    !             WRITE(numout,*) "update : SUM(fmmap(",ip,jp,")) = ",sum_fm
    !          ENDIF
    !       ENDIF
    !    ENDDO
    ! ENDDO
    !
    !
    ! The number of maximum vegetation map points in the GCM grid should
    ! also be computed and not imposed here.
    !
    nbvmax = 200
    !
    callsign="Spinup clearcut map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(npts, nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in sub_index','','')
       ENDIF
       sub_index(:,:,:)=0

       ALLOC_ERR=-1
       ALLOCATE(sub_area(npts, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_spinup_clearcut', 'error in sub_area','','') 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(npts, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_ful, lat_ful, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       !
       nbvmax = nbvmax * 2
    ENDDO

    ! Obtain the majority clearcut flag value after regrdding.
    DO ipts = 1, npts
       ! We need to do this for every PFT
       DO ivma=1,nvmap
          ! For this point, we need to see what dominant FM types are
          ! nearby.
          sumf=zero
          ! If it's not a forest, we don't care.
          IF(.NOT. is_tree(start_index(ivma)))THEN
             spinup_clearcut(ipts,ivma)=0
          ELSE
             DO ibvm=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ipts,ibvm) <= zero ) EXIT
                ip = sub_index(ipts,ibvm,1)
                jp = sub_index(ipts,ibvm,2)
                ifm=fmmap_i(ip,jp,ivma)
                sumf = sumf+ifm
                !WRITE(numout,*) 'jp: ',jp,'ip: ',ip,'ifm: ',ifm,'sumf: ',sumf,'ibmv: ',ibvm
             ENDDO
             ! As long as any sub-pixel goes through clearcut, the coarse pixel
             ! goes through clearcut.
             IF (sumf .GT. zero) THEN
                spinup_clearcut(ipts,ivma) = flag_spinup_clearcut
             ELSE
                spinup_clearcut(ipts,ivma) = zero
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    WRITE(numout,*) 'Done with interpolating the spinup clearcut map'

    !
    DEALLOCATE(fmmap_i)
    DEALLOCATE(fmmap_r)
    DEALLOCATE(lat_lu,lon_lu)
    DEALLOCATE(lat_ful,lon_ful)
    DEALLOCATE(mask)
    IF(ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF(ALLOCATED(sub_area)) DEALLOCATE(sub_area)

   IF ( printlev_loc >= 5 ) WRITE(numout,*) 'Leaving sapiens_forestry_read_spinup_clearcut'

  END SUBROUTINE sapiens_forestry_read_spinup_clearcut

!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_read_litter
!!
!>\BRIEF        Read in a map that gives the litter demand for each pixel.
!!
!! DESCRIPTION   :  Reads in a map for the litter demand of each pixel and interpolates
!!                  to find the demand for all pixels that we are interested in.  The map
!!                  in this case has units of gC/year.  We convert it to be gC/m**2/year, since
!!                  all the other units in the model are given per square meter.  The litter
!!                  demand maps can be calculated based on, for example, the demand needed
!!                  per head of livestock for storing them indoors in the winter.  This practice
!!                  reached its peak during the 19th century and rapidly declined afterwards
!!                  in Europe.  If you are not interested in historical simulations or improved 
!!                  estimation of soil carbon pools during a historical spinup, you probably 
!!                  don't need this.
!!     
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::litter_demand
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_read_litter ( npts, lalo, neighbours, resolution, contfrac, litter_demand )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    REAL(r_std), DIMENSION(:), INTENT(out)              :: litter_demand   !! The litter removed from each pixel due
                                                                           !! to litter raking at the end of the year.
                                                                           !! @tex $(gC year^{-1})$ @endtex

    !! 0.3 Modified fields


    !! 0.4 Local variables
    !
    ! PARAMETERS...taken from grid.f90.  I don't know why this variable is not in constants.f90.
    ! default resolution (m)
    REAL(r_std), PARAMETER :: default_resolution = 250000.

    CHARACTER(LEN=80)                                      :: filename        !! A string to hold the file name
    LOGICAL                                                :: debug=.FALSE.   !! A flag to print out debugging information.
    INTEGER(i_std)                                         :: fid             !! The ID of the NetCDF file.
    INTEGER(i_std)                                         :: nb_coord        !! The number of coordinates in the NetCDF file    
    INTEGER(i_std)                                         :: nb_gat          !!
    INTEGER(i_std)                                         :: nb_var          !! The number of variables in the NetCDF file
    INTEGER(i_std)                                         :: nb_dim          !! The number of dimensions in the NetCDF file
    INTEGER(i_std)                                         :: iml, jml, lml,ivmi,xxx !! indices
    INTEGER(i_std)                                         :: ip, inbv, jp    !! indices
    LOGICAL                                                :: l_ex            !! A flag which indicates if a variable
                                                                              !! exists in the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                 :: lat_lu, lon_lu  !! The latitude and longitude read in from
                                                                              !! the NetCDF file    
    INTEGER,DIMENSION(flio_max_var_dims)                   :: l_d_w           !! List of the dimension lengths of the variable
                                                                              !! in the NetCDF file
    INTEGER(i_std)                                         :: ipts            !! index
    INTEGER(i_std)                                         :: ALLOC_ERR       !! A flag tripped if we have an error in allocation
    REAL(r_std),DIMENSION(:,:),ALLOCATABLE                 :: littermap       !! The map read in from the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_ful, lon_ful
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: ibvm, nbvmax,ifm
    LOGICAL ::           ok_interpol ! optionnal return of aggregate_2d
    REAL(r_std)                                            :: area_sum, coslat
    REAL(r_std),DIMENSION(nfm_types)                       :: fm_sum
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: temp_resolution
    REAL(r_std)                :: fillvalue
    INTEGER :: i_rc,i_v,no_fill
    LOGICAL ::           ltemp ! temporary logical variable
!_ ================================================================================================================================

    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Entering sapiens_forestry_read_litter'

    litter_demand(:)=zero


   !
    !Config Key   = LITTER_FILE
    !Config Desc  = Name of file from which the litter raking map is to be read
    !Config If    = OK_STOMATE
    !Config Def   = litter_map.nc
    !Config Help  = 
    !Config Units = [FILE]
    !
    filename = 'litter_map.nc'
    CALL getin_p('LITTER_FILE',filename)

    WRITE(numout,*)'reading the litter_map: ',filename

    IF (is_root_prc) THEN
       IF (debug) THEN
          WRITE(numout,*) "Entering sapiens_forestry_read_litter. Debug mode."
          WRITE (*,'(/," --> fliodmpf")')
          CALL fliodmpf (TRIM(filename))
          WRITE (*,'(/," --> flioopfd")')
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (debug) THEN
          WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF
    CALL bcast(nb_coord)
    CALL bcast(nb_var)
    CALL bcast(nb_gat)

    ! I need the fill value of the map.  This doesn't seem to available through
    ! the IOISPL routines that I can see.
    ! NF90_INQ_VAR_FILL seems to be only valid for NetCDF-'?  So I am going to make
    ! a HUGE assumption here about the fill value of the maps.  I hope that none of the values
    ! in the maps are above this value.
    fillvalue=1e20
    WRITE(numout,*) "WARNING: You are using a litter demand map.  I cannot determine"
    WRITE(numout,*) "         the fill value of this map.  I am assuming that the fill"
    WRITE(numout,*) "         value is 1e20 and that all your real values are below that."
    WRITE(numout,*) "         CONFIRM THIS WITH YOUR LITTER MAP!"
    

    ! This finds the number of longitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    iml=l_d_w(1)
    WRITE(numout,*) "Litter map: iml =",iml

    ! This finds the number of latitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    jml=l_d_w(1)
    WRITE(numout,*) "Litter map: jml =",jml

    ! Now find the number of PFTs in the file.  If this is not equal to the number
    ! of PFTs that we actually have, that's a problem and we'll crash.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="LITTER_DEMAND",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)

    !

    ! Allocate the map that will be read in
    WRITE(numout,*) 'Reading the litter file',l_d_w
    !
    ALLOC_ERR=-1
    ALLOCATE(littermap(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of littermap : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in littermap','','')
    ENDIF
    ALLOC_ERR=-1

    ! This reads in the map that is in the file
    IF (is_root_prc) THEN
       CALL fliogetv (fid,"LITTER_DEMAND", littermap, start=(/ 1, 1 /), count=(/ iml, jml /))
    ENDIF

    CALL bcast(littermap)

    ! Now I need to get the latitude and longitude
    ! First, get the axes from the map file.
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in lat_lu','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in lon_lu','','')
    ENDIF
    IF (is_root_prc) THEN
       CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
    ENDIF
    ! Do we need to close the file??
    !ASLA IF (is_root_prc) THEN
    !ASLA   CALL flinclo(fid)
    !ASLA ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)

    ! Now I can interpolate.  This should be straightforward, since
    ! we are interpolating real numbers.
    
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in lat_ful','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in lon_ful','','')
    ENDIF
    !
    DO ip=1,iml
       lon_ful(ip,:)=lon_lu(ip)
    ENDDO
    DO jp=1,jml
       lat_ful(:,jp)=lat_lu(jp)
    ENDDO
    
    ! The map is in absolute numbers.  We cannot interpolate absolute numbers, we want
    ! to add them together (or subtract them, if we are going to finer resolution).  Therefore,
    ! let's convert all the map values to per m**2 first, interpolate that, and then
    ! convert back to absolute numbers.
    ! This is ugly.  Unfortunately, we had to steal this from grid.f90, since the resolution is
    ! not passed for every grid square, just those that we are interested in.
    ALLOC_ERR=-1
    ALLOCATE(temp_resolution(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of temp_resolution : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in temp_resolution','','')
    ENDIF
    DO ip=1,iml
       DO jp=1,jml
          !
          ! 2 resolution
          !
          
          !
          ! 2.1 longitude
          !
          
          ! prevent infinite resolution at the pole
          coslat = MAX( COS( lat_ful(ip,jp) * pi/180. ), mincos )     
          IF ( iml .GT. 1 ) THEN
             
             IF ( ip .EQ. 1 ) THEN
                temp_resolution(ip,jp,1) = &
                     ABS( lon_ful(ip+1,jp) - lon_ful(ip,jp) ) * &
                     pi/180. * R_Earth * coslat
             ELSEIF ( ip .EQ. iml ) THEN
                temp_resolution(ip,jp,1) = &
                     ABS( lon_ful(ip,jp) - lon_ful(ip-1,jp) ) * &
                     pi/180. * R_Earth * coslat
             ELSE
                temp_resolution(ip,jp,1) = &
                     ABS( lon_ful(ip+1,jp) - lon_ful(ip-1,jp) )/2. *&
                     pi/180. * R_Earth * coslat
             ENDIF
             
          ELSE
             
             temp_resolution(ip,jp,1) = default_resolution
             
          ENDIF

          !
          ! 2.2 latitude
          !
          
          IF ( jml .GT. 1 ) THEN

             IF ( jp .EQ. 1 ) THEN
                temp_resolution(ip,jp,2) = &
                     ABS( lat_ful(ip,jp) - lat_ful(ip,jp+1) ) * &
                     pi/180. * R_Earth
             ELSEIF ( jp .EQ. jml ) THEN
                temp_resolution(ip,jp,2) = &
                     ABS( lat_ful(ip,jp-1) - lat_ful(ip,jp) ) * &
                     pi/180. * R_Earth
             ELSE
                temp_resolution(ip,jp,2) = &
                     ABS( lat_ful(ip,jp-1) - lat_ful(ip,jp+1) )/2. *&
                     pi/180. * R_Earth
             ENDIF
             
          ELSE
             
             temp_resolution(ip,jp,2) = default_resolution

          ENDIF
             
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(littermap(ip,jp))
#else
                ltemp=isnan(littermap(ip,jp))
#endif
          IF( .NOT. ltemp)THEN

             IF(littermap(ip,jp) .LT. fillvalue)THEN
                ! Convert to per square meter for the interpolation if this is a real map value.
                littermap(ip,jp)=littermap(ip,jp)/temp_resolution(ip,jp,1)/temp_resolution(ip,jp,2)
             ENDIF
             
          ENDIF
          
       ENDDO

    ENDDO

    !
    WRITE(numout,*) 'Reading the LITTER file'
    !
    !

    !
    ! Mask of permitted variables.
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in mask','','')
    ENDIF
    !
    mask(:,:) = 0
    DO ip=1,iml
       DO jp=1,jml
          ! If the grid has any litter, we should interpolate
          ! this grid point.
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(littermap(ip,jp))
#else
                ltemp=isnan(littermap(ip,jp))
#endif
          IF( .NOT. ltemp)THEN
             IF ( littermap(ip,jp) .GT. min_sechiba .AND. littermap(ip,jp) .LT. fillvalue ) THEN
                mask(ip,jp) = 1
                IF (debug) THEN
                   WRITE(numout,*) "update : littermap(",ip,jp,")) = ",littermap(ip,jp)
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !
    ! The number of maximum vegetation map points in the GCM grid should
    ! also be computed and not imposed here.
    !
    nbvmax = 200
    !
    callsign="Litter map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(npts, nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in sub_index','','')
       ENDIF
       sub_index(:,:,:)=0

       ALLOC_ERR=-1
       ALLOCATE(sub_area(npts, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_litter', 'error in sub_area','','')
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(npts, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_ful, lat_ful, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       !
       nbvmax = nbvmax * 2
    ENDDO
    !
    DO ipts = 1, npts

       litter_demand(ipts)=zero
       area_sum=zero

       DO ibvm=1, nbvmax
          ! Leave the do loop if all sub areas are treated, sub_area <= 0
          IF ( sub_area(ipts,ibvm) <= zero ) EXIT
          ip = sub_index(ipts,ibvm,1)
          jp = sub_index(ipts,ibvm,2)
          area_sum=area_sum+sub_area(ipts,ibvm)
          litter_demand(ipts)=litter_demand(ipts)+sub_area(ipts,ibvm)*littermap(ip,jp)
       ENDDO

       IF(area_sum .LT. min_stomate)THEN
          WRITE(numout,*) 'Missing data for the litter map on this land point!'
          WRITE(numout,*) 'ipts,lalo',ipts,lalo(ipts,:)
          litter_demand(ipts)=zero
       ELSE
          litter_demand(ipts)=litter_demand(ipts)/area_sum
          ! This is in gC/year/m**2 from the interpolation.  I need to convert it back
          ! to absolute numbers over the ORCHIDEE cell (area)
          litter_demand(ipts)=litter_demand(ipts)*area(ipts)
       ENDIF

    ENDDO

    !
    DEALLOCATE(littermap)
    DEALLOCATE(temp_resolution)
    DEALLOCATE(lat_lu,lon_lu)
    DEALLOCATE(lat_ful,lon_ful)
    DEALLOCATE(mask)
    IF(ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF(ALLOCATED(sub_area)) DEALLOCATE(sub_area)


   IF ( printlev_loc >= 5 ) WRITE(numout,*) 'Leaving sapiens_forestry_read_litter'

  END SUBROUTINE sapiens_forestry_read_litter

!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_litter_raking
!!
!>\BRIEF       Redistribute litter between PFTs due to litter raking 
!!
!! DESCRIPTION   :  
!!                  (1) The transfer of forest litter to crop lands within each grid cell depends on the litter demands. 
!!                  This is detemined by carbon litter demand maps read in stomate.f90.
!!                  (2) Determine the CN ratio of litter pools, because the CN
!!                  ratio should be the same before and after the litter raking.
!!                  (3) Determine the amount of litter removed from each PFT. We already know 
!!                  the total amount af litter to be transfered from the forest 
!!                  to the crops, but it has to be distibuted amongst the PTF. This is determinde 
!!                  by the fraction of above ground litter of a given PFT to the total above ground 
!!                  litter in the grid cell. (The litter fraction of a PFT can very well be 
!!                  different from veget_max). 
!!                  (4) Remove nitrogen from the forest litter pools, such that a
!!                  constant CN ratio is kept before and after the litter raking.
!!                  (5) Move the litter of C and N to the crop PFTs in the grid. 
!!     
!!                  Be obs on the units in this subroutine - the changes between
!!                  mass per grid cell and mass per m2, because the litter
!!                  demand maps gC per grid cell. 
!!
!! RECENT CHANGE(S) : Nitrogen was added to the litter raking subroutine during
!!                    summer 2019
!!
!! MAIN OUTPUT VARIABLE(S): ::litter pools are modified
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_litter_raking ( npts, veget_max, resolution, litter, lrake_frac )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: npts            !! Domain size - number of pixels 
                                                                          !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: veget_max       !! "Maximal" coverage fraction of a PFT (LAI 
                                                                          !! -> infinity) on ground
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: resolution      !! The size in m of each grid-box in X and Y

    !! 0.2 Output
    REAL(r_std), DIMENSION(:,:), INTENT(out)           :: lrake_frac      !! Relative amount of litter that raked (-)

    !! 0.3 Modified fields
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter          !! Metabolic and structural litter, above 
                                                                          !! and below ground 
                                                                          !! @tex $(gC m^{-2})$ @endtex 

    !! 0.4 Local variables
    INTEGER(i_std)                                     :: ipts, ivm       !! Indices
    REAL(r_std)                                        :: total_litter    !! Total amount of aboveground forest litter 
                                                                          !! on the pixel (gC)
    REAL(r_std)                                        :: ind_litter      !! Resulting litter per PFT
    REAL(r_std)                                        :: crop_frac       !! Fraction of croplands on the pixel (-)
    REAL(r_std)                                        :: forest_frac     !! Fraction of forests on the pixel (-)
    REAL(r_std)                                        :: litter_transfer !! Litter that we would like to move (demand)
                                                                          !! from the forests to the crops
                                                                          !! @tex $(gC year^{-1})$ @endtex
    REAL(r_std),DIMENSION(nlitt)                       :: litter_ratio    !! Litter in a PFT as a fraction of all the 
                                                                          !! forest litter in the pixel (-)
    REAL(r_std),DIMENSION(nvm)                         :: lfrac           !! ??? how does this differ from litter_ratio
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern    !! Contains the components of the internal
                                                                          !! mass balance chech for this routine
                                                                          !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: closure_intern  !! Check closure of internal mass balance
                                                                          !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_start      !! Start pool of this routine 
                                                                          !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_end        !! End pool of this routine 
                                                                          !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std)                                        :: temp_sum        !! Temporary storage
    INTEGER(i_std)                                     :: iele, ilitt     !! Indices
    INTEGER(i_std)                                     :: ilevs, imbc     !! Indices
    REAL(r_std), DIMENSION(npts,nvm)                   :: pft_area        !! The absolute area covered by this PFT
                                                                          !! @tex $(m^2)$ @endtex
    REAL(r_std)                                        :: litter_lost     !! Litter removed from the forest PFTs
    REAL(r_std)                                        :: litter_gained   !! Litter added to the crop PFTs
    LOGICAL                                            :: lall_litter     !! Flag to catch numerical issues
    REAL(r_std), DIMENSION(npts,nvm)                   :: veget_max_begin !! temporary storage of veget_max to check area conservation
    REAL(r_std), DIMENSION(nlitt,nvm)                  :: CN_litter       !! CN ratios for the above ground litter pools 
                                                                          !! before  litter raking 
    REAL(r_std), DIMENSION(nlitt,nvm,nelements)        :: litter_lost_pft !! The litter lost due to litter raking for each PFT (gC or gN per grid cell)
    REAL(r_std), DIMENSION(npts)                       :: litter_transfer_N !! The total tranfers of N from the forested PFT to the crops (gN per grid cell).
    REAL(r_std), DIMENSION(npts)                       :: litter_demand   !! The litter removed from each pixel due to litter raking at the end of the year (gC year^{-1})
    REAL(r_std)                                        :: ind_litter_N    !! Resulting N litter per PFT (gN per m2)
!_ ================================================================================================================================

    IF (firstcall_sapiens_forestry) THEN
       ! Initialize local printlev if it is not already done
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_forestry=.FALSE.
    END IF

    IF ( printlev_loc >= 3 ) WRITE(numout,*) 'Entering sapiens_forestry_litter_raking'

    
    !! 1.1 Read fluxes from litter raking from file
    !      Note, this subroutine is only called once a year, that's why it can be done in the main loop
    CALL sapiens_forestry_read_litter(npts, lalo, neighbours, resolution, &
         contfrac, litter_demand)
 
    !! 1.2 Initialize check for mass balance closure
    
    !  Initial values of litter
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
    
       DO iele = 1,nelements
       
          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilitt = 1,nlitt
             DO ilevs = 1,nlevs
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     litter(:,ilitt,:,ilevs,iele) * veget_max(:,:)
             ENDDO
          ENDDO
       
       ENDDO ! # nelements

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1

    !! 1.4 Initialize
    lrake_frac(:,:) = zero
    litter_transfer_N(:) = zero

    !! Surface areas in m2
    DO ipts=1,npts
       pft_area(ipts,:)=area(ipts)* veget_max(ipts,:)
    ENDDO

    ! We need to find how much litter we can get from the forests for each 
    ! grid point.
    DO ipts=1,npts

       litter_lost=zero
       litter_gained=zero
       total_litter=zero
       CN_litter(:,:) = zero
       litter_lost_pft(:,:,:) = zero
     

       DO ivm=1,nvm
          IF(is_tree(ivm))THEN
             ! We will take all the above ground litter, since that is what 
             ! a farmer could pick up.  Unlikely they dug into the soil to 
             ! collect litter. Notice these are in absolute units, not per 
             ! square meter.
             total_litter=total_litter+SUM(litter(ipts,:,ivm,iabove,icarbon))&
                  *pft_area(ipts,ivm)
          ENDIF
       ENDDO
 
       ! We need to determine the CN ratio of the above grpund litter pools
       ! such that we keep the CN ratios constant before and after the litter
       ! raking
       DO ivm = 1,nvm
          IF(is_tree(ivm)) THEN

             IF(SUM(litter(ipts,:,ivm,iabove,icarbon)) .GT. min_stomate &
                .AND. SUM(litter(ipts,:,ivm,iabove,initrogen)) .GT. min_stomate)THEN 
            
               ! In case of empty litter N by component   
               WHERE(litter(ipts,:,ivm,iabove,initrogen) .GT. zero)   
                     CN_litter(:,ivm) =litter(ipts,:,ivm,iabove,icarbon)/ &
                              &litter(ipts,:,ivm,iabove,initrogen)
               ENDWHERE

               IF(ivm .EQ. test_pft .AND. printlev_loc.GE.4) THEN
                  WRITE(numout,*)'ivm',ivm
                  WRITE(numout,*)'CN_litter before',CN_litter(:,ivm)
               ENDIF

             ELSEIF(SUM(litter(ipts,:,ivm,iabove,icarbon)) .LT. min_stomate &
                   .AND. SUM(litter(ipts,:,ivm,iabove,initrogen)) .LT.min_stomate)THEN

             ! If there us no litter, we cannot take the CN ratio 
           
             ELSE

             ! There is something wrong with the CN ratio of the litter pools
                WRITE(numout,*)'ivm',ivm
                WRITE(numout,*)'litter C',litter(ipts,:,ivm,iabove,icarbon)
                WRITE(numout,*)'litter N',litter(ipts,:,ivm,iabove,initrogen)   
                CALL ipslerr_p (3,'ERROR: in litter raking','Something is wrong with the CN ratio of litter &
                             pools', 'One pools is empty, while the other is not','')
             ENDIF
          ENDIF
       ENDDO

       ! If we have no forest litter, we cannot rake anything.
       IF(total_litter .LT. min_stomate)CYCLE

       ! The forests are tricky. The litter they produce will not 
       ! necessarily be proportional to the area they take up.  
       ! We want to reduce the amount of litter in each forest in a 
       ! way that doesn't leave any litter pools negative.
       lfrac(:)=zero
       DO ivm=1,nvm
          IF(is_tree(ivm))THEN
             lfrac(ivm)=SUM(litter(ipts,:,ivm,iabove,icarbon))&
                  *pft_area(ipts,ivm)/total_litter
          ENDIF
       ENDDO

       ! Now compute the amount of the litter that we need to take to 
       ! satisfy our demand.
       litter_transfer=litter_demand(ipts)

       ! We don't take more litter than we have. Putting in 
       ! a logical flag to prevent us from having numerical 
       ! issues.
       lall_litter=.FALSE.
       IF(litter_transfer .GT. total_litter) THEN
          lall_litter=.TRUE.
          litter_transfer=total_litter
       ENDIF

       ! How many crop PFTs are present on this grid square?
       crop_frac=zero
       DO ivm=1,nvm
          IF(.NOT. natural(ivm))THEN
             crop_frac=crop_frac+veget_max(ipts,ivm)
          ENDIF
       ENDDO

       ! Litter raking is only conducted if we have crops in the grid cell.
       IF(crop_frac .GT. min_stomate)THEN
          DO ivm=1,nvm

             IF(veget_max(ipts,ivm) .LE. min_stomate) CYCLE

             ! If we have a forest, we remove litter.  Keep the ratio of 
             ! the litter pools the same.
             IF(is_tree(ivm))THEN
                ind_litter=SUM(litter(ipts,:,ivm,iabove,icarbon))*&
                     pft_area(ipts,ivm)

                ! I guess there is a chance that a forest has no litter.  
                ! This means there is nothing to take away.
                IF(ind_litter .LE. min_stomate)CYCLE

                litter_ratio(:)=litter(ipts,:,ivm,iabove,icarbon)*&
                     pft_area(ipts,ivm)/ind_litter
                lrake_frac(ipts,ivm)=lfrac(ivm)*litter_transfer/ind_litter

                ! We don't want to have negative litter values
                IF(lall_litter)THEN
                   ind_litter=zero
                ELSE
                   ind_litter=ind_litter-lfrac(ivm)*litter_transfer
                ENDIF
                litter_lost=litter_lost-lfrac(ivm)*litter_transfer

                ! We need to estimate the lost litter of C for each PFT to
                ! calculate lost N as well. But we also need to estimate, to
                ! total N lost. The below is total estimates per grid cell

                litter_lost_pft(:,ivm,icarbon)=-lfrac(ivm)*litter_transfer*&
                                               litter_ratio(:)
                WHERE(CN_litter(:,ivm) .GT. zero)
                    litter_lost_pft(:,ivm,initrogen)=-lfrac(ivm)*litter_transfer*&
                                               litter_ratio(:)/CN_litter(:,ivm)
                ENDWHERE
                litter_transfer_N(ipts)=litter_transfer_N(ipts)-&
                                        SUM(litter_lost_pft(:,ivm,initrogen))  

                ! This is a problem.  We are taking differences of numbers
                ! which are 1e9 or so, which means we could have a negative 
                ! value that is insignificant. I put the flag to try to 
                ! prevent that.
                IF(.NOT. lall_litter .AND. ind_litter .LT. -min_stomate)THEN
                   WRITE(numout,*) "Took away too much litter!"
                   WRITE(numout,*) "ipts,ivm:",ipts,ivm
                   WRITE(numout,*) "lfrac(ivm),ind_litter,litter_transfer: ",&
                        lfrac(ivm),SUM(litter(ipts,:,ivm,iabove,icarbon))*&
                        pft_area(ipts,ivm),litter_transfer
                ENDIF
                litter(ipts,:,ivm,iabove,icarbon)=litter_ratio(:)*ind_litter&
                     /pft_area(ipts,ivm)
                WHERE(CN_litter(:,ivm) .GT. zero)
                    litter(ipts,:,ivm,iabove,initrogen)=litter_ratio(:)*ind_litter&
                         /pft_area(ipts,ivm)/CN_litter(:,ivm) 
                ENDWHERE
                ! The resulting litter pool can also be determined be
                ! subtraction now that the litter lost per pft has been
                ! calculated (litter_lost_pft) to keep constant CN ratio. It
                ! might be more intuitive when people are looking into the code.
                ! It has been verifyed that both methods gives the same results.
                ! litter(ipts,:,ivm,iabove,initrogen)=litter(ipts,:,ivm,iabove,initrogen)+&
                !           litter_lost_pft(:,ivm,initrogen)/pft_area(ipts,ivm)

                !--- DEBUG ---! 
                IF (printlev_loc>=4) THEN
                   IF(ipts .EQ. test_grid .and.ivm .EQ. test_pft) THEN
                      WRITE(numout,*)'ivm',ivm
                      WRITE(numout,*)'litter',litter(ipts,:,ivm,iabove,initrogen)
                      WRITE(numout,*)'CNratio after',litter(ipts,:,ivm,iabove,icarbon)/&
                                           litter(ipts,:,ivm,iabove,initrogen)
                      WRITE(numout,*)'litter_lost_pft carbon',litter_lost_pft(:,ivm,icarbon)
                      WRITE(numout,*)'litter_lost_pft nitrogen',litter_lost_pft(:,ivm,initrogen)
                      WRITE(numout,*)'litter_transfer', litter_transfer
                      WRITE(numout,*)'lfrac', lfrac(ivm)
                      WRITE(numout,*)'litter_ratio',litter_ratio
                   ENDIF 
                ENDIF  

             ELSEIF( .NOT. natural(ivm))THEN
 
                ! If we have a crop, we add litter.
                ind_litter=SUM(litter(ipts,:,ivm,iabove,icarbon))* &
                     pft_area(ipts,ivm)
                ind_litter_N=SUM(litter(ipts,:,ivm,iabove,initrogen))* &
                     pft_area(ipts,ivm)  

                ! Is it possible that we have no crop litter but we 
                ! have veget_max?  Maybe.
                IF(ind_litter .LT. min_stomate)THEN
                   litter_ratio(:)=un/REAL(nlitt)
                ELSE
                   litter_ratio(:)=litter(ipts,:,ivm,iabove,icarbon)*&
                        pft_area(ipts,ivm)/ind_litter
                ENDIF
                ! This is the total litter we have for this PFT in 
                ! this pixel
                ind_litter=ind_litter+veget_max(ipts,ivm)/crop_frac*&
                     litter_transfer
                ind_litter_N=ind_litter_N+veget_max(ipts,ivm)/crop_frac*&
                             litter_transfer_N(ipts)
                ! could add a dimension to ind_litter, in stead of making a new
                ! variable.
                litter_gained=litter_gained+veget_max(ipts,ivm)/&
                     crop_frac*litter_transfer
                litter(ipts,:,ivm,iabove,icarbon)=litter_ratio(:)*&
                     ind_litter/pft_area(ipts,ivm)
                litter(ipts,:,ivm,iabove,initrogen)=litter_ratio(:)*&
                     ind_litter_N/pft_area(ipts,ivm)

 
             ENDIF ! checking to see if this PFT is a crop or forest

          ENDDO ! loop over PFTs

       ELSE

          ! We have no agricultural PFTs to move the litter to, so we
          ! assume that no litter raking was done in this grid square.

       ENDIF ! checking to see if we have crops

    ENDDO ! loop over pixels
    
 !! 2. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 2.2 Check surface area
       CALL check_vegetation_area("stomate_prescribe", npts, veget_max_begin, &
            veget_max,'pft')

       ! 2.3 Mass balance closure 
       ! 2.3.1 Calculate final biomass
       pool_end = zero
       DO iele = 1,nelements
          DO ilitt = 1,nlitt
             DO ilevs = 1,nlevs
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     litter(:,ilitt,:,ilevs,iele) * veget_max(:,:)
             ENDDO
          ENDDO
       ENDDO

       !! 2.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,ipoolchange,iele) = -un * &
               (pool_end(:,:,iele) - pool_start(:,:,iele))
       ENDDO

       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(:,test_pft,imbc,iele)
             !-
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("sapiens_forestry_litter_raking", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'pixel')

    ENDIF ! err_act.GT.1
    
    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Leaving sapiens_forestry_litter_raking'
    
  END SUBROUTINE sapiens_forestry_litter_raking

  
!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_set_species_change
!!
!>\BRIEF         Set species_change_map that gives the new PFT that needs to be planted
!!               when the current PFT is harvested.
!!               This is done by reading a map or a parameter from run.def
!!
!! DESCRIPTION   :  
!!
!!  If we are doing species change, we need to read in the map or force
!!  the PFT values.
!!  When we implement the species change we usually want to change FM as
!!  well (and sometime we even need it, i.e.,
!!  when the current FM strategy is to coppice and we change the species
!!  to a conifer tree than we will need to change FM because conifers are
!!  typically not coppiced). Species change and FM change were implemented
!!  separatly because that was easier to test and debug but not all 
!!  combinations were tested. The recommended setting is read de FM_desired
!!  map when using species change.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): species_change_map
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_set_species_change ( npts, lalo, neighbours, resolution, contfrac, species_change_map )

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(out)    :: species_change_map !! The number of a PFT that each PFT will
                                                                              !! be converted into after death.
    !! 0.3 Modified fields
    !! 0.4 Local variables
    INTEGER(i_std)                                      :: ivm
    INTEGER(i_std), DIMENSION(npts,nvmap)               :: species_map_temp   !! Temporary file


    ! Initialize printlev_loc
    IF (firstcall_sapiens_forestry) THEN
       ! Initialize local printlev if it is not already done
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_forestry=.FALSE.
    END IF
    
    ! Read species change map
    IF (ok_read_species_change_map)THEN
       ! Read a map
       IF (printlev_loc >= 3) WRITE(numout,*) 'In sapiens_forestry_set_species_change, set species_change_map by reading a map'
       
       CALL sapiens_forestry_read_species_change(npts, lalo, neighbours, &
            resolution, contfrac, species_map_temp)
       
       ! If we are using age classes, we read in the map in the same
       ! way but then we change it a bit to account for age classes.
       IF (nagec .GT. 1)THEN
          ! All age classes of the same PFT will have the same management
          DO ivm = 1,nvm
             species_change_map(:,ivm)=species_map_temp(:,agec_group(ivm))
          ENDDO
       ELSE
          ! The number of pfts on the map is identical
          ! to the number of PFTs used in the simulation
          WRITE(numout,*) 'Reading single value for species change'
          species_change_map(:,:)=species_map_temp(:,:)
       END IF
    ELSE
       ! Read values from run.def
       IF (printlev_loc >= 3) WRITE(numout,*) 'In sapiens_forestry_set_species_change, set species_change_map by reading values from run.def'
       IF (printlev_loc >= 3) WRITE(numout,*) 'Use species_change_force= ', species_change_force
       IF (species_change_force .EQ. -9999) THEN
          ! If we end-up here the user has set-up the model
          ! such that we use the species_change code but only
          ! to change the management. To keep the code simple
          ! we will use a dummy species_change_map
          
          ! All age classes of the same PFT will have the same
          ! management
          DO ivm = 1,nvm
             species_change_map(:,ivm) = agec_group(ivm)
          ENDDO
       ELSE
          ! Use a single value (::species_change_force) instead 
          ! of the information from a map. This was implemented
          ! as a feature for testing and/or debugging. It allows
          ! to test on a single pixel without reading maps.
          species_change_map(:,:)=species_change_force
       END IF
    ENDIF ! end ok_read_species_change_map
    
          
  END SUBROUTINE sapiens_forestry_set_species_change

  
!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_read_species_change
!!
!>\BRIEF        Read in a map that gives the PFT that needs to be planted
!!               when the current s PFT is harvested.
!!
!! DESCRIPTION   :  
!!
!!       NOTE: This routine was mostly copied from slowproc where the PFTmap is read in.
!!             Grid interpolation is used, but only to look at the nearby pixels to see
!!             see which management strategy is dominant.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::species_map
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_read_species_change ( npts, lalo, neighbours, resolution, contfrac, species_map )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)         :: species_map     !! The number of a PFT that each PFT will
                                                                           !! be converted into after death.

    !! 0.3 Modified fields


    !! 0.4 Local variables
    CHARACTER(LEN=80)                                      :: filename        !! A string to hold the file name
    LOGICAL                                                :: debug=.FALSE.   !! A flag to print out debugging information.
    INTEGER(i_std)                                         :: fid             !! The ID of the NetCDF file.
    INTEGER(i_std)                                         :: nb_coord        !! The number of coordinates in the NetCDF file    
    INTEGER(i_std)                                         :: nb_gat          !!
    INTEGER(i_std)                                         :: nb_var          !! The number of variables in the NetCDF file
    INTEGER(i_std)                                         :: nb_dim          !! The number of dimensions in the NetCDF file
    INTEGER(i_std)                                         :: iml, jml, lml,ivm !! indices
    INTEGER(i_std)                                         :: ip, inbv, jp    !! indices
    LOGICAL                                                :: l_ex            !! A flag which indicates if a variable
                                                                              !! exists in the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                 :: lat_lu, lon_lu  !! The latitude and longitude read in from
                                                                              !! the NetCDF file    
    INTEGER,DIMENSION(flio_max_var_dims)                   :: l_d_w           !! List of the dimension lengths of the variable
                                                                              !! in the NetCDF file
    INTEGER(i_std)                                         :: ipts            !! index
    INTEGER(i_std)                                         :: ALLOC_ERR       !! A flag tripped if we have an error in allocation
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE               :: scmap_r         !! The map read in from the NetCDF file
    INTEGER(i_std),DIMENSION(:,:,:),ALLOCATABLE            :: scmap_i         !! The integer form of the map read in
    INTEGER(i_std)                                         :: closest_lat     !! The index of the closest latitude we found.
    INTEGER(i_std)                                         :: closest_lon     !! The index of the closest longitude we found.
    REAL(r_std)                                            :: distance        !! The distance from the current point to the 
                                                                              !! point on the map
    REAL(r_std)                                            :: closest_dist    !! The distance to the closet point we've found.
    INTEGER                                                :: large_int       !! A number which indicates that the grid data
                                                                              !! is not available
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_ful, lon_ful
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: ibvm, nbvmax,ifm, ivma
    LOGICAL ::           ok_interpol ! optionnal return of aggregate_2d
    REAL(r_std)                                            :: sum_sc,sumf
    REAL(r_std),DIMENSION(nvmap)                       :: sc_sum
    LOGICAL ::           ltemp ! temporary logical variable
  !_ ================================================================================================================================


    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Entering sapiens_forestry_read_species_change'

    ! This integer has to be large enough that it never shows up on the map, without being
    ! so large that is causes overflows.  This needs to be larger than the largest possible
    ! value multiplied by the number of PFTs.  Since each value should never be larger than
    ! the number of PFTs, this should suffice.
    large_int = nvmap*nvmap+1

   !
    !Config Key   = SPECIES_CHANGE_FILE
    !Config Desc  = Name of file from which the species change map is to be read
    !Config If    = OK_STOMATE
    !Config Def   = replant_species.nc
    !Config Help  = The name of the file to be opened to read the species to replant
    !Config         for each pixel and PFT.
    !Config Units = [FILE]
    !
    filename = 'replant_species.nc'
    CALL getin_p('SPECIES_CHANGE_FILE',filename)

    IF (is_root_prc) THEN
       IF (debug) THEN
          WRITE(numout,*) "Entering sapiens_forestry_read_species_change. Debug mode."
          WRITE (*,'(/," --> fliodmpf")')
          CALL fliodmpf (TRIM(filename))
          WRITE (*,'(/," --> flioopfd")')
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (debug) THEN
          WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF
    CALL bcast(nb_coord)
    CALL bcast(nb_var)
    CALL bcast(nb_gat)

    ! This finds the number of longitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    iml=l_d_w(1)
    WRITE(numout,*) "Species change Map: iml =",iml

    ! This finds the number of latitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    jml=l_d_w(1)
    WRITE(numout,*) "Species change Map: jml =",jml

    ! Now find the number of PFTs in the file.  If this is not equal to the number
    ! of PFTs that we actually have, that's a problem and we'll crash.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="NEWSPECIES",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    lml=l_d_w(3)

    IF (lml /= nvmap) THEN
       WRITE(numout,*) 'lml = ',lml
       WRITE(numout,*) 'nvmap = ',nvmap
       WRITE(numout,*) 'Stopping. '
       CALL ipslerr_p (3,'sapians_forestry', &
            &          'Problem with species change map.','lml /= nvmap', &
            &          '(number of pft must be equal)')
    ENDIF
    !

    ! Allocate the map that will be read in
    WRITE(numout,*) 'Reading the species change map'
    !
    ALLOC_ERR=-1
    ALLOCATE(scmap_r(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of scmap_r : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in scmap_r','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(scmap_i(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of scmap_i : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in scmap_i','','') 
    ENDIF

    ! This reads in the map that is in the file
    IF (is_root_prc) THEN
       scmap_r(:,:,:)=large_int*2.0
       CALL fliogetv (fid,"NEWSPECIES", scmap_r, start=(/ 1, 1, 1 /), count=(/ iml, jml, nvmap /))
       ! Right now the values are all real, but they should be integers.
       ! Careful, NINT might not work if the precision of scmap_r is not single.
       ! In that case, IDNINT should be used.
       DO ip=1,iml
          DO jp=1,jml
             DO ivma=1,nvmap
                !  There will be some fill values in here.  If we pass
                ! a large fill value to NINT, it crashes.  So let's 
                ! test for it
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(scmap_r(ip,jp,ivma))    
#else
                ltemp=isnan(scmap_r(ip,jp,ivma))    
#endif
                IF(ltemp)THEN
                   scmap_i(ip,jp,ivma)=large_int

                ELSE
                   IF(scmap_r(ip,jp,ivma) .GE. 0.0 .AND. scmap_r(ip,jp,ivma) < large_int)THEN
                      scmap_i(ip,jp,ivma)=NINT(scmap_r(ip,jp,ivma))
                   ELSE
                      ! This value should be big enough that we don't barl ourselves
                      ! below.
                      scmap_i(ip,jp,ivma)=large_int
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL bcast(scmap_i)

    ! Now I need to get the latitude and longitude
    ! First, get the axes from the map file.
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in lat_lu','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in lon_lu','','')
    ENDIF
    IF (is_root_prc) THEN
       CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)

    ! Now I can interpolate.  Remember that we are dealing with integer
    ! values and PFTs.  Therefore, we cannot do a strict
    ! interpolation, since that might gives values of the PFT to be
    ! 2.3, or something ridiculous.  We cannot just round the number, either,
    ! since we could have some squares with poplar (PFT 15) and oak (PFT 13), 
    ! rounding would give a completely different PFT (14).  So we look to see what the 
    ! most prevalent type of PFT of nearby pixels is, and then just take that.
    
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in lat_ful','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in lon_ful','','')
    ENDIF
    !
    DO ip=1,iml
       lon_ful(ip,:)=lon_lu(ip)
    ENDDO
    DO jp=1,jml
       lat_ful(:,jp)=lat_lu(jp)
    ENDDO
    !
    !
    WRITE(numout,*) 'Reading the SPECIES CHANGE file'
    !
    !

    !
    ! Mask of permitted variables.
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in mask','','')
    ENDIF
    !
    mask(:,:) = 0
    DO ip=1,iml
       DO jp=1,jml
          ! If at least one PFT has a management strategy, we should interpolate
          ! this grid point.
          sum_sc=SUM(scmap_i(ip,jp,:))
          IF ( sum_sc .GT. min_sechiba .AND. sum_sc .LT. large_int) THEN
             mask(ip,jp) = 1
             IF (debug) THEN
                WRITE(numout,*) "update : SUM(scmap(",ip,jp,")) = ",sum_sc
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !
    ! The number of maximum vegetation map points in the GCM grid should
    ! also be computed and not imposed here.
    !
    nbvmax = 200
    !
    callsign="Species change map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(npts, nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in sub_index','','')
       ENDIF
       sub_index(:,:,:)=0

       ALLOC_ERR=-1
       ALLOCATE(sub_area(npts, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_species_change', 'error in sub_area','','') 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(npts, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_ful, lat_ful, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       !
       nbvmax = nbvmax * 2
    ENDDO
    !
    DO ipts = 1, npts
       ! For this point, we need to see what dominant new species are
       ! nearby.
       sumf=zero
       ! We need to do this for every PFT
       DO ivma=1,nvmap
          ! If it's not a forest, we don't care.  We will
          ! always replant it as the same PFT.
          IF(.NOT. is_tree(start_index(ivma)))THEN
             species_map(ipts,ivma)=ivma
          ELSE
             sc_sum(:)=0.0
             DO ibvm=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ipts,ibvm) <= zero ) EXIT
                ip = sub_index(ipts,ibvm,1)
                jp = sub_index(ipts,ibvm,2)
                ifm=scmap_i(ip,jp,ivma)
                sc_sum(ifm)=sc_sum(ifm)+sub_area(ipts,ibvm)
             ENDDO
             ! Whichever FM type has the most area, we use that one.
             species_map(ipts,ivma)=MAXLOC(sc_sum(:),1)
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE(scmap_i)
    DEALLOCATE(scmap_r)
    DEALLOCATE(lat_lu,lon_lu)
    DEALLOCATE(lat_ful,lon_ful)
    DEALLOCATE(mask)
    IF(ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF(ALLOCATED(sub_area)) DEALLOCATE(sub_area)


    IF ( printlev_loc >= 5 ) WRITE(numout,*) 'Leaving sapiens_forestry_read_species_change'

  END SUBROUTINE sapiens_forestry_read_species_change


!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_set_desired_fm
!!
!>\BRIEF        Read in a map that gives the desired forest management strategy
!!              after a new PFT was planted.
!!
!! DESCRIPTION  
!!     
!!
!! FM = 1 : No human intervention (ORCHIDEE default)
!!      2 : Thinnings based on the RDI, clearcuts based on tree density,
!!              annual increment, and tree diameter.  Thinnings from above and
!!              from below are determined by the sign of thstrat.
!!      3 : Coppices
!!      4 : Short rotation coppices
!!      5 : Continuous cover forestry
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::forest_managed
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_set_desired_fm ( npts, lalo, neighbours, resolution, &
       contfrac, forest_managed, fm_change_map )

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
    !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)     :: forest_managed  !! Forest management flag
    
    !! 0.2 Output
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(out)    :: fm_change_map   !! A map which gives the desired FM strategy when
                                                                           !! the PFT will be replanted after a clearcut.

    !! 0.3 Modified fields
    !! 0.4 Local variables
    INTEGER(i_std)                                      :: ivm
    INTEGER(i_std), DIMENSION(npts,nvmap)               :: fm_map_temp     !! Temporary file

    ! Initialize printlev_loc
    IF (firstcall_sapiens_forestry) THEN
       ! Initialize local printlev if it is not already done
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_forestry=.FALSE.
    END IF
       
    ! Read species change map
    IF (ok_read_desired_fm_map)THEN
       IF (printlev_loc>=4) WRITE(numout,*) 'In sapiens_forestry_set_desired_fm, set fm_change_map by reading a map'         
       ! If we are using age classes, we read in the map in the same
       ! way but then we change it a bit to account for age classes.
       CALL sapiens_forestry_read_desired_fm(npts, lalo, neighbours, &
               resolution, contfrac, fm_map_temp)
       
       IF (nagec .GT. 1)THEN
          ! All age classes of the same PFT will have the same management
          DO ivm = 1,nvm
             fm_change_map(:,ivm)=fm_map_temp(:,agec_group(ivm))
          ENDDO
       ELSE
          ! The number of pfts on the map is identical
          ! to the number of PFTs used in the simulation
          fm_change_map(:,:)=fm_map_temp(:,:)
       ENDIF
    ELSE
       IF (printlev_loc>=4) WRITE(numout,*) 'In sapiens_forestry_set_desired_fm, set fm_change_map with values from run.def'
       IF (printlev_loc>=4) WRITE(numout,*) 'Use fm_change_force= ',fm_change_force
       IF (fm_change_force .EQ. -9999) THEN
          ! If we end-up here the user has set-up the model
          ! such that we use the species_change code but only
          ! to change the species while keeping the management
          ! constant. To keep the code simple we will use a 
          ! dummy fm_change_map
          fm_change_map(:,:)=forest_managed(:,:)
       ELSE
          ! Use a single value (::species_change_force) instead 
          ! of the information from a map. This was implemented
          ! as a feature for testing and/or debugging. It allows
          ! to test on a single pixel without reading maps.
          fm_change_map(:,:)=fm_change_force
       ENDIF
    ENDIF
    
    
  END SUBROUTINE sapiens_forestry_set_desired_fm

    
!! ================================================================================================================================
!! SUBROUTINE    : sapiens_forestry_read_desired_fm
!!
!>\BRIEF        Read in a map that gives the desired forest management strategy
!!              after a new PFT was planted.
!!
!! DESCRIPTION   :  
!!     
!!
!! FM = 1 : No human intervention (ORCHIDEE default)
!!      2 : Thinnings based on the RDI, clearcuts based on tree density,
!!              annual increment, and tree diameter.  Thinnings from above and
!!              from below are determined by the sign of thstrat.
!!      3 : Coppices
!!      4 : Short rotation coppices
!!      5 : Continuous cover forestry
!!
!!       NOTE: This routine was mostly copied from slowproc where the PFTmap is read in.
!!             Grid interpolation is used, but only to look at the nearby pixels to see
!!             see which management strategy is dominant.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::forest_managed
!!
!! REFERENCE(S)   : 
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE sapiens_forestry_read_desired_fm ( npts, lalo, neighbours, resolution, &
       contfrac, desired_managed )

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of pixels 
                                                                           !! (dimensionless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)          :: neighbours      !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: resolution      !! The size in m of each grid-box in X and Y
    REAL(r_std), DIMENSION(:), INTENT(in)               :: contfrac        !! Fraction of continent in the grid

    !! 0.2 Output
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)         :: desired_managed !! Forest management flag: 1 = unmanaged (orchidee
                                                                           !! standard), 2= self-thinning only, 3= coppices,
                                                                           !! 4= short rotation coppices, 5= continous
                                                                           !! cover forestry
    !! 0.3 Modified fields


    !! 0.4 Local variables
    CHARACTER(LEN=80)                                      :: filename        !! A string to hold the file name
    LOGICAL                                                :: debug=.FALSE.   !! A flag to print out debugging information.
    INTEGER(i_std)                                         :: fid             !! The ID of the NetCDF file.
    INTEGER(i_std)                                         :: nb_coord        !! The number of coordinates in the NetCDF file    
    INTEGER(i_std)                                         :: nb_gat          !!
    INTEGER(i_std)                                         :: nb_var          !! The number of variables in the NetCDF file
    INTEGER(i_std)                                         :: nb_dim          !! The number of dimensions in the NetCDF file
    INTEGER(i_std)                                         :: iml,jml,lml,ivm !! indices
    INTEGER(i_std)                                         :: ip, inbv, jp    !! indices
    LOGICAL                                                :: l_ex            !! A flag which indicates if a variable
                                                                              !! exists in the NetCDF file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)                 :: lat_lu, lon_lu  !! The latitude and longitude read in from
                                                                              !! the NetCDF file    
    INTEGER,DIMENSION(flio_max_var_dims)                   :: l_d_w           !! List of the dimension lengths of the variable
                                                                              !! in the NetCDF file
    INTEGER(i_std)                                         :: ipts            !! index
    INTEGER(i_std)                                         :: ALLOC_ERR       !! A flag tripped if we have an error in allocation
    REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE               :: fmmap_r         !! The map read in from the NetCDF file
    INTEGER(i_std),DIMENSION(:,:,:),ALLOCATABLE            :: fmmap_i         !! The integer form of the map read in
    INTEGER(i_std)                                         :: closest_lat     !! The index of the closest latitude we found.
    INTEGER(i_std)                                         :: closest_lon     !! The index of the closest longitude we found.
    REAL(r_std)                                            :: distance        !! The distance from the current point to the 
                                                                              !! point on the map
    REAL(r_std)                                            :: closest_dist    !! The distance to the closet point we've found.
    INTEGER                                                :: large_int       !! A number which indicates that the grid data
                                                                              !! is not available
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)               :: lat_ful, lon_ful
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)            :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)               :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)          :: sub_index
    CHARACTER(LEN=30)                                      :: callsign
    INTEGER(i_std)                                         :: ibvm, nbvmax,ifm, ivma
    LOGICAL                                                :: ok_interpol     !! optionnal return of aggregate_2d
    REAL(r_std)                                            :: sum_fm,sumf
    REAL(r_std),DIMENSION(nfm_types)                       :: fm_sum
    LOGICAL                                                :: ltemp           !! temporary logical variable
  !_ ================================================================================================================================

    IF ( printlev_loc >= 4 ) WRITE(numout,*) 'Entering sapiens_forestry_read_desired_fm'

    ! This integer has to be large enough that it never shows up on the map, without being
    ! so large that is causes overflows.  Since none of the points on the map should be
    ! larger than the number of FM strategies we have (nfm_types), this is sufficiently big.
    large_int = nvmap*nfm_types+1

    !Config Key   = FM_DESIRED_FILE
    !Config Desc  = Name of file from which the forest management map is to be read
    !Config If    = OK_STOMATE
    !Config Def   = FM_desired.nc
    !Config Help  = The name of the file to be opened to read a forest management
    !Config         map (including a layer for every PFT) is given here. 
    !Config Units = [FILE]
    filename = 'FM_desired.nc'
    CALL getin_p('FM_DESIRED_FILE',filename)

    IF (is_root_prc) THEN
       IF (debug) THEN
          WRITE(numout,*) "Entering sapiens_forestry_read_desired_fm. Debug mode."
          WRITE (*,'(/," --> fliodmpf")')
          CALL fliodmpf (TRIM(filename))
          WRITE (*,'(/," --> flioopfd")')
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (debug) THEN
          WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF
    CALL bcast(nb_coord)
    CALL bcast(nb_var)
    CALL bcast(nb_gat)

    ! This finds the number of longitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    iml=l_d_w(1)
    WRITE(numout,*) "FM desired Map: iml =",iml

    ! This finds the number of latitude points in the file.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    jml=l_d_w(1)
    WRITE(numout,*) "FM desired Map: jml =",jml

    ! Now find the number of PFTs in the file.  If this is not equal to the number
    ! of PFTs that we actually have, that's a problem and we'll crash.
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="FM_STRAT",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    lml=l_d_w(3)

    IF (lml /= nvmap) THEN
       WRITE(numout,*) 'lml = ',lml
       WRITE(numout,*) 'nvmap = ',nvmap
       WRITE(numout,*) 'Stopping. '
       CALL ipslerr_p (3,'sapiens_forestry', &
            &          'Problem with the desired forest management strateg map.',&
            &          'lml /= nvmap', '(number of pft must be equal)')
    ENDIF

    ! Allocate the map that will be read in
    WRITE(numout,*) 'Reading the Desired Forest Management strategy file'
    ALLOC_ERR=-1
    ALLOCATE(fmmap_r(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of fmmap_r : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in fmmap_r','','') 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(fmmap_i(iml,jml,nvmap), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of fmmap_i : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in fmmpap_i','','')
    ENDIF

    ! This reads in the map that is in the file
    IF (is_root_prc) THEN
       fmmap_r(:,:,:)=large_int*2.0
       CALL fliogetv (fid,"FM_STRAT", fmmap_r, start=(/ 1, 1, 1 /), &
            count=(/ iml, jml, nvmap /))

       ! Right now the values are all real, but they should be integers.
       ! Careful, NINT might not work if the precision of fmmap_r is not 
       ! single. In that case, IDNINT should be used.
       DO ip=1,iml
          DO jp=1,jml
             DO ivma=1,nvmap

                ! There will be some fill values in here.  If we pass
                ! a large fill value to NINT, it crashes.  So let's 
                ! test for it
#ifdef CPP_IEEE_ARITHMETIC
                ltemp=IEEE_IS_NAN(fmmap_r(ip,jp,ivma))    
#else
                ltemp=isnan(fmmap_r(ip,jp,ivma))    
#endif
                IF(ltemp)THEN
                   fmmap_i(ip,jp,ivma)=large_int
                ELSE
                   IF(fmmap_r(ip,jp,ivma) .GE. 0.0 .AND. fmmap_r(ip,jp,ivma) < large_int) THEN
                      fmmap_i(ip,jp,ivma) = NINT(fmmap_r(ip,jp,ivma))
                   ELSE

                      ! This value should be big enough that we don't barl ourselves
                      ! below.
                      fmmap_i(ip,jp,ivma)=large_int
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    CALL bcast(fmmap_i)

    ! Now I need to get the latitude and longitude
    ! First, get the axes from the map file.
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in lat_lu','','') 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in lon_lu','','') 
    ENDIF
    IF (is_root_prc) THEN
       CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)

    ! Now I can interpolate.  Remember that we are dealing with integer
    ! values and management strategies.  Therefore, we cannot do a strict
    ! interpolation, since that might gives values of the strategy to be
    ! 2.3, or something ridiculous.  We cannot just round the number, either,
    ! since we could have some squares with 1 (no management) and some with
    ! 3 (coppices);  An average of that would be 2 (high stands), which doesn't
    ! make any sense.  So we look to see what the most prevalent type of
    ! management of nearby pixels is, and then just take that.
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in lat_ful','','')
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in lon_ful','','')
    ENDIF
    
    DO ip=1,iml
       lon_ful(ip,:)=lon_lu(ip)
    ENDDO
    DO jp=1,jml
       lat_ful(:,jp)=lat_lu(jp)
    ENDDO
    
    ! Mask of permitted variables.
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in mask','','')
    ENDIF
    
    mask(:,:) = 0
    DO ip=1,iml
       DO jp=1,jml
          ! If at least one PFT has a management strategy, we should interpolate
          ! this grid point.
          sum_fm=SUM(fmmap_i(ip,jp,:))
          IF ( sum_fm .GT. min_sechiba .AND. sum_fm .LT. large_int) THEN
             mask(ip,jp) = 1
             IF (debug) THEN
                WRITE(numout,*) "update : SUM(fmmap(",ip,jp,")) = ",sum_fm
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    
    
    ! The number of maximum vegetation map points in the GCM grid should
    ! also be computed and not imposed here.
    nbvmax = 200
    
    callsign="Desired Forest Management map"
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       
       ALLOC_ERR=-1
       ALLOCATE(sub_index(npts, nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in sub_index','','')
       ENDIF
       sub_index(:,:,:)=0

       ALLOC_ERR=-1
       ALLOCATE(sub_area(npts, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          CALL ipslerr_p (3,'sapiens_forestry_read_desired_fm', 'error in sub_area','','') 
       ENDIF
       sub_area(:,:)=zero
       
       CALL aggregate_p(npts, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_ful, lat_ful, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       
       nbvmax = nbvmax * 2
    ENDDO
    
    DO ipts = 1, npts

       ! For this point, we need to see what dominant FM types are
       ! nearby.
       sumf=zero

       ! We need to do this for every PFT
       DO ivma=1,nvmap

          IF(.NOT. is_tree(start_index(ivma)))THEN

             ! If it's not a forest, we don't care.
             desired_managed(ipts,ivma)=0

          ELSE

             fm_sum(:)=0.0
             DO ibvm=1, nbvmax
                ! Leave the do loop if all sub areas are treated, 
                ! sub_area <= 0
                IF ( sub_area(ipts,ibvm) <= zero ) EXIT
                ip = sub_index(ipts,ibvm,1)
                jp = sub_index(ipts,ibvm,2)
                ifm=fmmap_i(ip,jp,ivma)
                fm_sum(ifm)=fm_sum(ifm)+sub_area(ipts,ibvm)
             ENDDO

             ! Whichever FM type has the most area, we use that one.
             desired_managed(ipts,ivma)=MAXLOC(fm_sum(:),1)
          ENDIF
       ENDDO
    ENDDO
    
    DEALLOCATE(fmmap_i)
    DEALLOCATE(fmmap_r)
    DEALLOCATE(lat_lu,lon_lu)
    DEALLOCATE(lat_ful,lon_ful)
    DEALLOCATE(mask)
    IF(ALLOCATED(sub_index)) DEALLOCATE(sub_index)
    IF(ALLOCATED(sub_area)) DEALLOCATE(sub_area)

    IF ( printlev_loc >= 5 ) WRITE(numout,*) 'Leaving sapiens_forestry_read_desired_fm'

  END SUBROUTINE sapiens_forestry_read_desired_fm

END MODULE sapiens_forestry
