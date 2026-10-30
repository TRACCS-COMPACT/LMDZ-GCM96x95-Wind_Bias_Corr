!
! $Id$
!
MODULE guide_loc_mod

!=======================================================================
!   Auteur:  F.Hourdin
!            F. Codron 01/09
!=======================================================================

  USE getparam, only: ini_getparam, fin_getparam, getpar
  USE Write_Field_loc
  USE netcdf, ONLY: nf90_nowrite, nf90_open, nf90_inq_varid, nf90_close, &
          nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_dimid, &
          nf90_inquire_dimension, nf90_enddef, nf90_def_dim, nf90_put_var, nf90_noerr, nf90_close, nf90_inq_varid, &
          nf90_redef, nf90_write, nf90_unlimited, nf90_float, nf90_clobber, nf90_64bit_offset, nf90_float, &
          nf90_create, nf90_def_var, nf90_open
  USE parallel_lmdz
  USE pres2lev_mod, only: pres2lev

  USE paramet_mod_h
IMPLICIT NONE

! ---------------------------------------------
! Declarations des cles logiques et parametres
! ---------------------------------------------
  INTEGER, PRIVATE, SAVE  :: iguide_read,iguide_int,iguide_sav
  INTEGER, PRIVATE, SAVE  :: nlevnc, guide_plevs
  LOGICAL, PRIVATE, SAVE  :: guide_u,guide_v,guide_T,guide_Q,guide_P
  LOGICAL, PRIVATE, SAVE  :: guide_hr,guide_teta
  LOGICAL, PRIVATE, SAVE  :: guide_BL,guide_reg,guide_add,gamma4,guide_zon
  LOGICAL, PRIVATE, SAVE  :: invert_p,invert_y,ini_anal
  LOGICAL, PRIVATE, SAVE  :: guide_2D,guide_sav,guide_modele
!FC
  LOGICAL, PRIVATE, SAVE  :: convert_Pa

  REAL, PRIVATE, SAVE     :: tau_min_u,tau_max_u
  REAL, PRIVATE, SAVE     :: tau_min_v,tau_max_v
  REAL, PRIVATE, SAVE     :: tau_min_T,tau_max_T
  REAL, PRIVATE, SAVE     :: tau_min_Q,tau_max_Q
  REAL, PRIVATE, SAVE     :: tau_min_P,tau_max_P

  REAL, PRIVATE, SAVE     :: lat_min_g,lat_max_g
  REAL, PRIVATE, SAVE     :: lon_min_g,lon_max_g
  REAL, PRIVATE, SAVE     :: tau_lon,tau_lat

  REAL, PRIVATE, SAVE     :: plim_guide_BL

  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: alpha_u,alpha_v
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: alpha_T,alpha_Q
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: alpha_P,alpha_pcor

! ---------------------------------------------
! Variables de guidage
! ---------------------------------------------
! Variables des fichiers de guidage
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: unat1,unat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: vnat1,vnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: tnat1,tnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: qnat1,qnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: pnat1,pnat2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: psnat1,psnat2
  REAL, ALLOCATABLE, DIMENSION(:),     PRIVATE, SAVE   :: apnc,bpnc
! Variables aux dimensions du modele
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: ugui1,ugui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: vgui1,vgui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: tgui1,tgui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: qgui1,qgui2
  REAL, ALLOCATABLE, DIMENSION(:),   PRIVATE, SAVE   :: psgui1,psgui2

  INTEGER,SAVE,PRIVATE :: ijbu,ijbv,ijeu,ijev !,ijnu,ijnv
  INTEGER,SAVE,PRIVATE :: jjbu,jjbv,jjeu,jjev,jjnu,jjnv


CONTAINS
!=======================================================================

  SUBROUTINE guide_init

    USE control_mod, ONLY: day_step
    USE serre_mod, ONLY: grossismx

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE pycpl, ONLY: lk_pycpl
IMPLICIT NONE




    INTEGER                :: error,ncidpl,rid,rcod
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_init'
    CHARACTER (len = 20)   :: namedim

! ---------------------------------------------
! Lecture des parametres:
! ---------------------------------------------
    call ini_getparam("nudging_parameters_out.txt")
! Variables guidees
    CALL getpar('guide_u',.true.,guide_u,'guidage de u')
    CALL getpar('guide_v',.true.,guide_v,'guidage de v')
    CALL getpar('guide_T',.true.,guide_T,'guidage de T')
    CALL getpar('guide_P',.true.,guide_P,'guidage de P')
    CALL getpar('guide_Q',.true.,guide_Q,'guidage de Q')
    CALL getpar('guide_hr',.true.,guide_hr,'guidage de Q par H.R')
    CALL getpar('guide_teta',.false.,guide_teta,'guidage de T par Teta')

    CALL getpar('guide_add',.false.,guide_add,'for�age constant?')
    CALL getpar('guide_zon',.false.,guide_zon,'guidage moy zonale')
    if (guide_zon .and. abs(grossismx - 1.) > 0.01) &
         call abort_gcm("guide_init", &
         "zonal nudging requires grid regular in longitude", 1)

!   Constantes de rappel. Unite : fraction de jour
    CALL getpar('tau_min_u',0.02,tau_min_u,'Cste de rappel min, u')
    CALL getpar('tau_max_u', 10.,tau_max_u,'Cste de rappel max, u')
    CALL getpar('tau_min_v',0.02,tau_min_v,'Cste de rappel min, v')
    CALL getpar('tau_max_v', 10.,tau_max_v,'Cste de rappel max, v')
    CALL getpar('tau_min_T',0.02,tau_min_T,'Cste de rappel min, T')
    CALL getpar('tau_max_T', 10.,tau_max_T,'Cste de rappel max, T')
    CALL getpar('tau_min_Q',0.02,tau_min_Q,'Cste de rappel min, Q')
    CALL getpar('tau_max_Q', 10.,tau_max_Q,'Cste de rappel max, Q')
    CALL getpar('tau_min_P',0.02,tau_min_P,'Cste de rappel min, P')
    CALL getpar('tau_max_P', 10.,tau_max_P,'Cste de rappel max, P')
    CALL getpar('gamma4',.false.,gamma4,'Zone sans rappel elargie')
    CALL getpar('guide_BL',.true.,guide_BL,'guidage dans C.Lim')
    CALL getpar('plim_guide_BL',85000.,plim_guide_BL,'BL top presnivs value')

! Sauvegarde du for�age
    CALL getpar('guide_sav',.false.,guide_sav,'sauvegarde guidage')
    CALL getpar('iguide_sav',4,iguide_sav,'freq. sauvegarde guidage')
    ! frequences f>0: fx/jour; f<0: tous les f jours; f=0: 1 seule fois.
    IF (iguide_sav.GT.0) THEN
       iguide_sav=day_step/iguide_sav
    ELSE if (iguide_sav == 0) then
       iguide_sav = huge(0)
    ELSE
       iguide_sav=day_step*iguide_sav
    ENDIF

! Guidage regional seulement (sinon constant ou suivant le zoom)
    CALL getpar('guide_reg',.false.,guide_reg,'guidage regional')
    CALL getpar('lat_min_g',-90.,lat_min_g,'Latitude mini guidage ')
    CALL getpar('lat_max_g', 90.,lat_max_g,'Latitude maxi guidage ')
    CALL getpar('lon_min_g',-180.,lon_min_g,'longitude mini guidage ')
    CALL getpar('lon_max_g', 180.,lon_max_g,'longitude maxi guidage ')
    CALL getpar('tau_lat', 5.,tau_lat,'raideur lat guide regional ')
    CALL getpar('tau_lon', 5.,tau_lon,'raideur lon guide regional ')

! Parametres pour lecture des fichiers
    CALL getpar('iguide_read',4,iguide_read,'freq. lecture guidage')
    CALL getpar('iguide_int',4,iguide_int,'freq. interpolation vert')
    IF (iguide_int.EQ.0) THEN
        iguide_int=1
    ELSEIF (iguide_int.GT.0) THEN
        iguide_int=day_step/iguide_int
    ELSE
        iguide_int=day_step*iguide_int
    ENDIF
    CALL getpar('guide_plevs',0,guide_plevs,'niveaux pression fichiers guidage')
    ! Pour compatibilite avec ancienne version avec guide_modele
    CALL getpar('guide_modele',.false.,guide_modele,'niveaux pression ap+bp*psol')
    IF (guide_modele) THEN
        guide_plevs=1
    ENDIF
!FC
    CALL getpar('convert_Pa',.true.,convert_Pa,'Convert Pressure levels in Pa')
    ! Fin raccord
    CALL getpar('ini_anal',.false.,ini_anal,'Etat initial = analyse')
    CALL getpar('guide_invertp',.true.,invert_p,'niveaux p inverses')
    CALL getpar('guide_inverty',.true.,invert_y,'inversion N-S')
    CALL getpar('guide_2D',.false.,guide_2D,'fichier guidage lat-P')

    call fin_getparam

! ---------------------------------------------
! Determination du nombre de niveaux verticaux
! des fichiers guidage
! ---------------------------------------------
    ncidpl=-99
    if (guide_plevs.EQ.1) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('apbp.nc',Nf90_NOWRITe, ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file apbp.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_plevs.EQ.2) then
       if (ncidpl.EQ.-99) then
          rcod=nf90_open('P.nc',Nf90_NOWRITe,ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file P.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif

    elseif (guide_u) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('u.nc',Nf90_NOWRITe,ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file u.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif

       endif


    elseif (guide_v) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('v.nc',nf90_nowrite,ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file v.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif


    elseif (guide_T) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('T.nc',nf90_nowrite,ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file T.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif



    elseif (guide_Q) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('hur.nc',nf90_nowrite, ncidpl)
          if (rcod.NE.nf90_noerr) THEN
             abort_message=' Nudging error -> no file hur.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif


    endif
    error=nf90_inq_dimid(ncidpl,'LEVEL',rid)
    IF (error.NE.nf90_noerr) error=nf90_inq_dimid(ncidpl,'PRESSURE',rid)
    IF (error.NE.nf90_noerr) THEN
        abort_message='Nudging: error reading pressure levels'
        CALL abort_gcm(modname,abort_message,1)
    ENDIF
    error=nf90_inquire_dimension(ncidpl,rid,len=nlevnc)
    write(*,*)trim(modname)//' : number of vertical levels nlevnc', nlevnc
    rcod = nf90_close(ncidpl)

    ! Python coupling : les champs python sont sur les niveaux du modele,
    ! l'ecrasement de unat2/vnat2 n'est possible que si nlevnc == llm
    IF (lk_pycpl .AND. nlevnc.NE.llm) THEN
        abort_message='pycpl: nlevnc /= llm, cannot overwrite unat2/vnat2'
        CALL abort_gcm(modname,abort_message,1)
    ENDIF

! ---------------------------------------------
! Allocation des variables
! ---------------------------------------------
    abort_message='nudging allocation error'

    ALLOCATE(apnc(nlevnc), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(bpnc(nlevnc), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    apnc=0.;bpnc=0.

    ALLOCATE(alpha_pcor(llm), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_u(ijb_u:ije_u), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_v(ijb_v:ije_v), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_T(ijb_u:ije_u), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_Q(ijb_u:ije_u), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_P(ijb_u:ije_u), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    alpha_u=0.;alpha_v=0;alpha_T=0;alpha_Q=0;alpha_P=0

    IF (guide_u) THEN
        ALLOCATE(unat1(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(ugui1(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(unat2(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(ugui2(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        unat1=0.;unat2=0.;ugui1=0.;ugui2=0.
    ENDIF

    IF (guide_T) THEN
        ALLOCATE(tnat1(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tgui1(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tnat2(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tgui2(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        tnat1=0.;tnat2=0.;tgui1=0.;tgui2=0.
    ENDIF

    IF (guide_Q) THEN
        ALLOCATE(qnat1(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qgui1(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qnat2(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qgui2(ijb_u:ije_u,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        qnat1=0.;qnat2=0.;qgui1=0.;qgui2=0.
    ENDIF

    IF (guide_v) THEN
        ALLOCATE(vnat1(iip1,jjb_v:jje_v,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vgui1(ijb_v:ije_v,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vnat2(iip1,jjb_v:jje_v,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vgui2(ijb_v:ije_v,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        vnat1=0.;vnat2=0.;vgui1=0.;vgui2=0.
    ENDIF

    IF (guide_plevs.EQ.2) THEN
        ALLOCATE(pnat1(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(pnat2(iip1,jjb_u:jje_u,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        pnat1=0.;pnat2=0.;
    ENDIF

    IF (guide_P.OR.guide_plevs.EQ.1) THEN
        ALLOCATE(psnat1(iip1,jjb_u:jje_u), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(psnat2(iip1,jjb_u:jje_u), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        psnat1=0.;psnat2=0.;
    ENDIF
    IF (guide_P) THEN
        ALLOCATE(psgui2(ijb_u:ije_u), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(psgui1(ijb_u:ije_u), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        psgui1=0.;psgui2=0.
    ENDIF

! ---------------------------------------------
!   Lecture du premier etat de guidage.
! ---------------------------------------------
    IF (guide_2D) THEN
        CALL guide_read2D(1)
    ELSE
        CALL guide_read(1)
    ENDIF
    IF (guide_v) vnat1=vnat2
    IF (guide_u) unat1=unat2
    IF (guide_T) tnat1=tnat2
    IF (guide_Q) qnat1=qnat2
    IF (guide_plevs.EQ.2) pnat1=pnat2
    IF (guide_P.OR.guide_plevs.EQ.1) psnat1=psnat2

  END SUBROUTINE guide_init

!=======================================================================
  SUBROUTINE guide_main(itau,ucov,vcov,teta,q,masse,ps)
    use exner_hyb_loc_m, only: exner_hyb_loc
    use exner_milieu_loc_m, only: exner_milieu_loc
    USE parallel_lmdz
    USE control_mod
    USE write_field_loc
    USE comconst_mod, ONLY: cpp, daysec, dtvr, kappa
    USE comvert_mod, ONLY: ap, bp, preff, presnivs, pressure_exner

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE




    ! Variables entree
    INTEGER,                           INTENT(IN)    :: itau !pas de temps
    REAL, DIMENSION (ijb_u:ije_u,llm), INTENT(INOUT) :: ucov,teta,q,masse
    REAL, DIMENSION (ijb_v:ije_v,llm), INTENT(INOUT) :: vcov
    REAL, DIMENSION (ijb_u:ije_u),     INTENT(INOUT) :: ps

    ! Variables locales
    LOGICAL, SAVE :: first=.TRUE.
!$OMP THREADPRIVATE(first)
    LOGICAL       :: f_out ! sortie guidage
    REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: f_addu ! var aux: champ de guidage
    REAL, ALLOCATABLE, SAVE, DIMENSION (:,:) :: f_addv ! var aux: champ de guidage
    ! Variables pour fonction Exner (P milieu couche)
    REAL, ALLOCATABLE, SAVE, DIMENSION (:,:,:)    :: pk
    REAL, ALLOCATABLE, SAVE, DIMENSION (:,:)        :: pks
    REAL                               :: unskap
    REAL, ALLOCATABLE, SAVE, DIMENSION (:,:)    :: p ! besoin si guide_P
    ! Compteurs temps:
    INTEGER, SAVE :: step_rea,count_no_rea,itau_test ! lecture guidage
!$OMP THREADPRIVATE(step_rea,count_no_rea,itau_test)
    REAL          :: ditau, dday_step
    REAL          :: tau,reste ! position entre 2 etats de guidage
    REAL, SAVE    :: factt ! pas de temps en fraction de jour
!$OMP THREADPRIVATE(factt)

    INTEGER       :: i,j,l
    CHARACTER(LEN=20) :: modname="guide_main"

!$OMP MASTER
    ijbu=ij_begin ; ijeu=ij_end
    jjbu=jj_begin ; jjeu=jj_end ; jjnu=jjeu-jjbu+1
    ijbv=ij_begin ; ijev=ij_end
    jjbv=jj_begin ; jjev=jj_end ; jjnv=jjev-jjbv+1
    IF (pole_sud) THEN
      ijeu=ij_end-iip1
      ijev=ij_end-iip1
      jjev=jj_end-1
      jjnv=jjev-jjbv+1
    ENDIF
    IF (pole_nord) THEN
      ijbu=ij_begin+iip1
      ijbv=ij_begin
    ENDIF
!$OMP END MASTER
!$OMP BARRIER

!    PRINT *,'---> on rentre dans guide_main'
!    CALL AllGather_Field(ucov,ip1jmp1,llm)
!    CALL AllGather_Field(vcov,ip1jm,llm)
!    CALL AllGather_Field(teta,ip1jmp1,llm)
!    CALL AllGather_Field(ps,ip1jmp1,1)
!    CALL AllGather_Field(q,ip1jmp1,llm)

!-----------------------------------------------------------------------
! Initialisations au premier passage
!-----------------------------------------------------------------------

    IF (first) THEN
        first=.FALSE.
!$OMP MASTER
        ALLOCATE(f_addu(ijb_u:ije_u,llm) )
        ALLOCATE(f_addv(ijb_v:ije_v,llm) )
        ALLOCATE(pk(iip1,jjb_u:jje_u,llm)  )
        ALLOCATE(pks(iip1,jjb_u:jje_u)  )
        ALLOCATE(p(ijb_u:ije_u,llmp1) )
        CALL guide_init
!$OMP END MASTER
!$OMP BARRIER
        itau_test=1001
        step_rea=1
        count_no_rea=0
! Calcul des constantes de rappel
        factt=dtvr*iperiod/daysec
!$OMP MASTER
        call tau2alpha(3, iip1, jjb_v, jje_v, factt, tau_min_v, tau_max_v, alpha_v)
        call tau2alpha(2, iip1, jjb_u, jje_u, factt, tau_min_u, tau_max_u, alpha_u)
        call tau2alpha(1, iip1, jjb_u, jje_u, factt, tau_min_T, tau_max_T, alpha_T)
        call tau2alpha(1, iip1, jjb_u, jje_u, factt, tau_min_P, tau_max_P, alpha_P)
        call tau2alpha(1, iip1, jjb_u, jje_u, factt, tau_min_Q, tau_max_Q, alpha_Q)
! correction de rappel dans couche limite
        if (guide_BL) then
             alpha_pcor(:)=1.
        else
            do l=1,llm
                alpha_pcor(l)=(1.+tanh(((plim_guide_BL-presnivs(l))/preff)/0.05))/2.
            enddo
        endif
!$OMP END MASTER
!$OMP BARRIER
! ini_anal: etat initial egal au guidage
        IF (ini_anal) THEN
            CALL guide_interp(ps,teta)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l=1,llm
              IF (guide_u) ucov(ijbu:ijeu,l)=ugui2(ijbu:ijeu,l)
              IF (guide_v) vcov(ijbv:ijev,l)=ugui2(ijbv:ijev,l)
              IF (guide_T) teta(ijbu:ijeu,l)=tgui2(ijbu:ijeu,l)
              IF (guide_Q) q(ijbu:ijeu,l)=qgui2(ijbu:ijeu,l)
            ENDDO

            IF (guide_P) THEN
!$OMP MASTER
                ps(ijbu:ijeu)=psgui2(ijbu:ijeu)
!$OMP END MASTER
!$OMP BARRIER
                CALL pression_loc(ijnb_u,ap,bp,ps,p)
                CALL massdair_loc(p,masse)
!$OMP BARRIER
            ENDIF
            RETURN
        ENDIF

    ENDIF !first

!-----------------------------------------------------------------------
! Lecture des fichiers de guidage ?
!-----------------------------------------------------------------------
    IF (iguide_read.NE.0) THEN
      ditau=real(itau)
      dday_step=real(day_step)
      IF (iguide_read.LT.0) THEN
          tau=ditau/dday_step/REAL(iguide_read)
      ELSE
          tau=REAL(iguide_read)*ditau/dday_step
      ENDIF
      reste=tau-AINT(tau)
      IF (reste.EQ.0.) THEN
          IF (itau_test.EQ.itau) THEN
            write(*,*)trim(modname)//' second pass in advreel at itau=',&
            itau
            CALL abort_gcm("guide_loc_lod","stopped",1)
          ELSE
!$OMP MASTER
              IF (guide_v) vnat1(:,jjbv:jjev,:)=vnat2(:,jjbv:jjev,:)
              IF (guide_u) unat1(:,jjbu:jjeu,:)=unat2(:,jjbu:jjeu,:)
              IF (guide_T) tnat1(:,jjbu:jjeu,:)=tnat2(:,jjbu:jjeu,:)
              IF (guide_Q) qnat1(:,jjbu:jjeu,:)=qnat2(:,jjbu:jjeu,:)
              IF (guide_plevs.EQ.2) pnat1(:,jjbu:jjeu,:)=pnat2(:,jjbu:jjeu,:)
              IF (guide_P.OR.guide_plevs.EQ.1) psnat1(:,jjbu:jjeu)=psnat2(:,jjbu:jjeu)
!$OMP END MASTER
!$OMP BARRIER
              step_rea=step_rea+1
              itau_test=itau
              if (is_master) then
                write(*,*)trim(modname)//' Reading nudging files, step ',&
                    step_rea,'after ',count_no_rea,' skips'
              endif
              IF (guide_2D) THEN
!$OMP MASTER
                  CALL guide_read2D(step_rea)
!$OMP END MASTER
!$OMP BARRIER
              ELSE
!$OMP MASTER
                  CALL guide_read(step_rea)
!$OMP END MASTER
!$OMP BARRIER
              ENDIF
              count_no_rea=0
          ENDIF
      ELSE
        count_no_rea=count_no_rea+1

      ENDIF
    ENDIF !iguide_read=0

!-----------------------------------------------------------------------
! Interpolation et conversion des champs de guidage
!-----------------------------------------------------------------------
    IF (MOD(itau,iguide_int).EQ.0) THEN
        CALL guide_interp(ps,teta)
    ENDIF
! Repartition entre 2 etats de guidage
    IF (iguide_read.NE.0) THEN
        tau=reste
    ELSE
        tau=1.
    ENDIF

!    CALL WriteField_u('ucov_guide',ucov)
!    CALL WriteField_v('vcov_guide',vcov)
!    CALL WriteField_u('teta_guide',teta)
!    CALL WriteField_u('masse_guide',masse)


!-----------------------------------------------------------------------
!   Ajout des champs de guidage
!-----------------------------------------------------------------------
! Sauvegarde du guidage?
    f_out=((MOD(itau,iguide_sav).EQ.0).AND.guide_sav)
    IF (f_out) THEN

!$OMP BARRIER
      CALL pression_loc(ijnb_u,ap,bp,ps,p)

!$OMP BARRIER
      if (pressure_exner) then
      CALL exner_hyb_loc( ijnb_u, ps, p, pks, pk)
      else
        CALL exner_milieu_loc( ijnb_u, ps, p, pks, pk )
      endif

!$OMP BARRIER

        unskap=1./kappa
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l = 1, llm
            DO j=jjbu,jjeu
                DO i =1, iip1
                    p(i+(j-1)*iip1,l) = preff * ( pk(i,j,l)/cpp) ** unskap
                ENDDO
            ENDDO
        ENDDO

        CALL guide_out("SP",jjp1,llm,p(ijb_u:ije_u,1:llm),1.)
    ENDIF

    if (guide_u) then
        if (guide_add) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
           f_addu(ijbu:ijeu,l)=(1.-tau)*ugui1(ijbu:ijeu,l)+tau*ugui2(ijbu:ijeu,l)
          ENDDO
        else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
           f_addu(ijbu:ijeu,l)=(1.-tau)*ugui1(ijbu:ijeu,l)+tau*ugui2(ijbu:ijeu,l)-ucov(ijbu:ijeu,l)
          ENDDO
        endif

!        CALL WriteField_u('f_addu',f_addu)

        if (guide_zon) CALL guide_zonave_u(1,llm,f_addu)
        CALL guide_addfield_u(llm,f_addu,alpha_u)
        IF (f_out) CALL guide_out("ua",jjp1,llm,(1.-tau)*ugui1(ijb_u:ije_u,:)+tau*ugui2(ijb_u:ije_u,:),factt)
        IF (f_out) CALL guide_out("u",jjp1,llm,ucov(ijb_u:ije_u,:),factt)
        IF (f_out) THEN
         ! Ehouarn: fill the gaps adequately...
         IF (ijbu>ijb_u) f_addu(ijb_u:ijbu-1,:)=0
         IF (ijeu<ije_u) f_addu(ijeu+1:ije_u,:)=0
         CALL guide_out("ucov",jjp1,llm,f_addu(ijb_u:ije_u,:)/factt,factt)
        ENDIF
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          ucov(ijbu:ijeu,l)=ucov(ijbu:ijeu,l)+f_addu(ijbu:ijeu,l)
        ENDDO

    endif

    if (guide_T) then
        if (guide_add) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
            f_addu(ijbu:ijeu,l)=(1.-tau)*tgui1(ijbu:ijeu,l)+tau*tgui2(ijbu:ijeu,l)
          ENDDO
        else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
           f_addu(ijbu:ijeu,l)=(1.-tau)*tgui1(ijbu:ijeu,l)+tau*tgui2(ijbu:ijeu,l)-teta(ijbu:ijeu,l)
          ENDDO
        endif
        if (guide_zon) CALL guide_zonave_u(2,llm,f_addu)
        CALL guide_addfield_u(llm,f_addu,alpha_T)
        IF (f_out) CALL guide_out("teta",jjp1,llm,f_addu(:,:)/factt,factt)
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          teta(ijbu:ijeu,l)=teta(ijbu:ijeu,l)+f_addu(ijbu:ijeu,l)
        ENDDO
    endif

    if (guide_P) then
        if (guide_add) then
!$OMP MASTER
            f_addu(ijbu:ijeu,1)=(1.-tau)*psgui1(ijbu:ijeu)+tau*psgui2(ijbu:ijeu)
!$OMP END MASTER
!$OMP BARRIER
        else
!$OMP MASTER
            f_addu(ijbu:ijeu,1)=(1.-tau)*psgui1(ijbu:ijeu)+tau*psgui2(ijbu:ijeu)-ps(ijbu:ijeu)
!$OMP END MASTER
!$OMP BARRIER
        endif
        if (guide_zon) CALL guide_zonave_u(2,1,f_addu(ijb_u:ije_u,1))
        CALL guide_addfield_u(1,f_addu(ijb_u:ije_u,1),alpha_P)
!       IF (f_out) CALL guide_out("ps",jjp1,1,f_addu(ijb_u:ije_u,1)/factt,factt)
!$OMP MASTER
        ps(ijbu:ijeu)=ps(ijbu:ijeu)+f_addu(ijbu:ijeu,1)
!$OMP END MASTER
!$OMP BARRIER
        CALL pression_loc(ijnb_u,ap,bp,ps,p)
        CALL massdair_loc(p,masse)
!$OMP BARRIER
    endif

    if (guide_Q) then
        if (guide_add) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
            f_addu(ijbu:ijeu,l)=(1.-tau)*qgui1(ijbu:ijeu,l)+tau*qgui2(ijbu:ijeu,l)
          ENDDO
        else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
            f_addu(ijbu:ijeu,l)=(1.-tau)*qgui1(ijbu:ijeu,l)+tau*qgui2(ijbu:ijeu,l)-q(ijbu:ijeu,l)
          ENDDO
        endif
        if (guide_zon) CALL guide_zonave_u(2,llm,f_addu)
        CALL guide_addfield_u(llm,f_addu,alpha_Q)
        IF (f_out) CALL guide_out("q",jjp1,llm,f_addu(:,:)/factt,factt)

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          q(ijbu:ijeu,l)=q(ijbu:ijeu,l)+f_addu(ijbu:ijeu,l)
        ENDDO
    endif

    if (guide_v) then
        if (guide_add) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
             f_addv(ijbv:ijev,l)=(1.-tau)*vgui1(ijbv:ijev,l)+tau*vgui2(ijbv:ijev,l)
          ENDDO

        else
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
            f_addv(ijbv:ijev,l)=(1.-tau)*vgui1(ijbv:ijev,l)+tau*vgui2(ijbv:ijev,l)-vcov(ijbv:ijev,l)
          ENDDO

        endif

        if (guide_zon) CALL guide_zonave_v(2,jjm,llm,f_addv(ijb_v:ije_v,:))

        CALL guide_addfield_v(llm,f_addv(ijb_v:ije_v,:),alpha_v)
        IF (f_out) CALL guide_out("v",jjm,llm,vcov(ijb_v:ije_v,:),factt)
        IF (f_out) CALL guide_out("va",jjm,llm,(1.-tau)*vgui1(ijb_v:ije_v,:)+tau*vgui2(ijb_v:ije_v,:),factt)
        IF (f_out) THEN
          ! Ehouarn: Fill in the gaps adequately
          IF (ijbv>ijb_v) f_addv(ijb_v:ijbv-1,:)=0
          IF (ijev<ije_v) f_addv(ijev+1:ije_v,:)=0
          CALL guide_out("vcov",jjm,llm,f_addv(ijb_v:ije_v,:)/factt,factt)
        ENDIF

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          vcov(ijbv:ijev,l)=vcov(ijbv:ijev,l)+f_addv(ijbv:ijev,l)
        ENDDO
    endif

  END SUBROUTINE guide_main


  SUBROUTINE guide_addfield_u(vsize,field,alpha)
! field1=a*field1+alpha*field2

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    ! input variables
    INTEGER,                      INTENT(IN)    :: vsize
    REAL, DIMENSION(ijb_u:ije_u),       INTENT(IN)    :: alpha
    REAL, DIMENSION(ijb_u:ije_u,vsize), INTENT(INOUT) :: field

    ! Local variables
    INTEGER :: l

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,vsize
      field(ijbu:ijeu,l)=alpha(ijbu:ijeu)*field(ijbu:ijeu,l)*alpha_pcor(l)
    ENDDO

  END SUBROUTINE guide_addfield_u


  SUBROUTINE guide_addfield_v(vsize,field,alpha)
! field1=a*field1+alpha*field2

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    ! input variables
    INTEGER,                      INTENT(IN)    :: vsize
    REAL, DIMENSION(ijb_v:ije_v),       INTENT(IN)    :: alpha
    REAL, DIMENSION(ijb_v:ije_v,vsize), INTENT(INOUT) :: field

    ! Local variables
    INTEGER :: l

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    DO l=1,vsize
      field(ijbv:ijev,l)=alpha(ijbv:ijev)*field(ijbv:ijev,l)*alpha_pcor(l)
    ENDDO

  END SUBROUTINE guide_addfield_v

!=======================================================================

  SUBROUTINE guide_zonave_u(typ,vsize,field)

    USE comgeom_mod_h
    USE comconst_mod, ONLY: pi

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
IMPLICIT NONE




    ! input/output variables
    INTEGER,                           INTENT(IN)    :: typ
    INTEGER,                           INTENT(IN)    :: vsize
    REAL, DIMENSION(ijb_u:ije_u,vsize), INTENT(INOUT) :: field

    ! Local variables
    LOGICAL, SAVE                :: first=.TRUE.
!$OMP THREADPRIVATE(first)

    INTEGER, DIMENSION (2), SAVE :: imin, imax ! averaging domain
!$OMP THREADPRIVATE(imin,imax)
    INTEGER                      :: i,j,l,ij
    REAL, DIMENSION (iip1)       :: lond       ! longitude in Deg.
    REAL, DIMENSION (jjb_u:jje_u,vsize):: fieldm     ! zon-averaged field

    IF (first) THEN
        first=.FALSE.
!Compute domain for averaging
        lond=rlonu*180./pi
        imin(1)=1;imax(1)=iip1;
        imin(2)=1;imax(2)=iip1;
        IF (guide_reg) THEN
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(1)=i
                IF (lond(i).LE.lon_max_g) imax(1)=i
            ENDDO
            lond=rlonv*180./pi
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(2)=i
                IF (lond(i).LE.lon_max_g) imax(2)=i
            ENDDO
        ENDIF
    ENDIF


!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,vsize
        fieldm(:,l)=0.
      ! Compute zonal average

!correction bug ici
! ---> a verifier
! ym         DO j=jjbv,jjev
         DO j=jjbu,jjeu
              DO i=imin(typ),imax(typ)
                  ij=(j-1)*iip1+i
                  fieldm(j,l)=fieldm(j,l)+field(ij,l)
              ENDDO
          ENDDO
          fieldm(:,l)=fieldm(:,l)/REAL(imax(typ)-imin(typ)+1)
    ! Compute forcing
          DO j=jjbu,jjeu
              DO i=1,iip1
                  ij=(j-1)*iip1+i
                  field(ij,l)=fieldm(j,l)
              ENDDO
          ENDDO
      ENDDO

  END SUBROUTINE guide_zonave_u


  SUBROUTINE guide_zonave_v(typ,hsize,vsize,field)

    USE comgeom_mod_h
    USE comconst_mod, ONLY: pi

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE




    ! input/output variables
    INTEGER,                           INTENT(IN)    :: typ
    INTEGER,                           INTENT(IN)    :: vsize
    INTEGER,                           INTENT(IN)    :: hsize
    REAL, DIMENSION(ijb_v:ije_v,vsize), INTENT(INOUT) :: field

    ! Local variables
    LOGICAL, SAVE                :: first=.TRUE.
!$OMP THREADPRIVATE(first)
    INTEGER, DIMENSION (2), SAVE :: imin, imax ! averaging domain
!$OMP THREADPRIVATE(imin, imax)
    INTEGER                      :: i,j,l,ij
    REAL, DIMENSION (iip1)       :: lond       ! longitude in Deg.
    REAL, DIMENSION (jjb_v:jjev,vsize):: fieldm     ! zon-averaged field

    IF (first) THEN
        first=.FALSE.
!Compute domain for averaging
        lond=rlonu*180./pi
        imin(1)=1;imax(1)=iip1;
        imin(2)=1;imax(2)=iip1;
        IF (guide_reg) THEN
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(1)=i
                IF (lond(i).LE.lon_max_g) imax(1)=i
            ENDDO
            lond=rlonv*180./pi
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(2)=i
                IF (lond(i).LE.lon_max_g) imax(2)=i
            ENDDO
        ENDIF
    ENDIF

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,vsize
      ! Compute zonal average
          fieldm(:,l)=0.
          DO j=jjbv,jjev
              DO i=imin(typ),imax(typ)
                  ij=(j-1)*iip1+i
                  fieldm(j,l)=fieldm(j,l)+field(ij,l)
              ENDDO
          ENDDO
          fieldm(:,l)=fieldm(:,l)/REAL(imax(typ)-imin(typ)+1)
    ! Compute forcing
          DO j=jjbv,jjev
              DO i=1,iip1
                  ij=(j-1)*iip1+i
                  field(ij,l)=fieldm(j,l)
              ENDDO
          ENDDO
      ENDDO


  END SUBROUTINE guide_zonave_v

!=======================================================================
  SUBROUTINE guide_interp(psi,teta)
    USE comgeom2_mod_h
  use exner_hyb_loc_m, only: exner_hyb_loc
    use exner_milieu_loc_m, only: exner_milieu_loc
  USE parallel_lmdz
  USE mod_hallo
  USE Bands
  USE comconst_mod, ONLY: cpp, kappa
  USE comvert_mod, ONLY: preff, pressure_exner, bp, ap, disvert_type
  USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE




  REAL, DIMENSION (iip1,jjb_u:jje_u),     INTENT(IN) :: psi ! Psol gcm
  REAL, DIMENSION (iip1,jjb_u:jje_u,llm), INTENT(IN) :: teta ! Temp. Pot. gcm

  LOGICAL, SAVE                      :: first=.TRUE.
!$OMP THREADPRIVATE(first)
  ! Variables pour niveaux pression:
  REAL, ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: plnc1,plnc2 !niveaux pression guidage
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)    :: plunc,plsnc !niveaux pression modele
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)     :: plvnc       !niveaux pression modele
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)  :: p           ! pression intercouches
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)    :: pls, pext   ! var intermediaire
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)    :: pbarx
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)     :: pbary
  ! Variables pour fonction Exner (P milieu couche)
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)    :: pk
  REAL ,ALLOCATABLE, SAVE, DIMENSION (:,:)        :: pks
  REAL                               :: unskap
  ! Pression de vapeur saturante
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:)      :: qsat
  !Variables intermediaires interpolation
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)    :: zu1,zu2
  REAL, ALLOCATABLE, SAVE,DIMENSION (:,:,:)     :: zv1,zv2

  INTEGER                            :: i,j,l,ij
  CHARACTER(LEN=20),PARAMETER :: modname="guide_interp"
  TYPE(Request),SAVE :: Req
!$OMP THREADPRIVATE(Req)

    if (is_master) write(*,*)trim(modname)//': interpolate nudging variables'
! -----------------------------------------------------------------
! Calcul des niveaux de pression champs guidage (pour T et Q)
! -----------------------------------------------------------------
    IF (first) THEN
!$OMP MASTER
      ALLOCATE(plnc1(iip1,jjb_u:jje_u,nlevnc) )
      ALLOCATE(plnc2(iip1,jjb_u:jje_u,nlevnc) )
      ALLOCATE(plunc(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(plsnc(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(plvnc(iip1,jjb_v:jje_v,llm) )
      ALLOCATE(p(iip1,jjb_u:jje_u,llmp1) )
      ALLOCATE(pls(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(pext(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(pbarx(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(pbary(iip1,jjb_v:jje_v,llm) )
      ALLOCATE(pk(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(pks (iip1,jjb_u:jje_u) )
      ALLOCATE(qsat(ijb_u:ije_u,llm) )
      ALLOCATE(zu1(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(zu2(iip1,jjb_u:jje_u,llm) )
      ALLOCATE(zv1(iip1,jjb_v:jje_v,llm) )
      ALLOCATE(zv2(iip1,jjb_v:jje_v,llm) )
!$OMP END MASTER
!$OMP BARRIER
    ENDIF




    IF (guide_plevs.EQ.0) THEN
!$OMP DO
        DO l=1,nlevnc
            DO j=jjbu,jjeu
                DO i=1,iip1
                    plnc2(i,j,l)=apnc(l)
                    plnc1(i,j,l)=apnc(l)
               ENDDO
            ENDDO
        ENDDO
    ENDIF

    if (first) then
        first=.FALSE.
!$OMP MASTER
        write(*,*)trim(modname)//' : check vertical level order'
        write(*,*)trim(modname)//' LMDZ :'
        do l=1,llm
          write(*,*)trim(modname)//' PL(',l,')=',(ap(l)+ap(l+1))/2. &
                  +psi(1,jjeu)*(bp(l)+bp(l+1))/2.
        enddo
        write(*,*)trim(modname)//' nudging file :'
        SELECT CASE (guide_plevs)
        CASE (0)
            do l=1,nlevnc
              write(*,*)trim(modname)//' PL(',l,')=',plnc2(1,jjbu,l)
            enddo
        CASE (1)
            DO l=1,nlevnc
              write(*,*)trim(modname)//' PL(',l,')=',&
                        apnc(l)+bpnc(l)*psnat2(1,jjbu)
            ENDDO
        CASE (2)
            do l=1,nlevnc
              write(*,*)trim(modname)//' PL(',l,')=',pnat2(1,jjbu,l)
            enddo
        END SELECT
        write(*,*)trim(modname)//' invert ordering: invert_p=',invert_p
        if (guide_u) then
            do l=1,nlevnc
              write(*,*)trim(modname)//' U(',l,')=',unat2(1,jjbu,l)
            enddo
        endif
        if (guide_T) then
            do l=1,nlevnc
              write(*,*)trim(modname)//' T(',l,')=',tnat2(1,jjbu,l)
            enddo
        endif
!$OMP END MASTER
    endif ! of if (first)

    if (guide_plevs /= 1 .or. guide_t .and. .not. guide_teta &
         .or. guide_q .and. guide_hr) then
       CALL pression_loc( ijnb_u, ap, bp, psi, p )
       if (disvert_type==1) then
          CALL exner_hyb_loc(ijnb_u,psi,p,pks,pk)
       else ! we assume that we are in the disvert_type==2 case
          CALL exner_milieu_loc(ijnb_u,psi,p,pks,pk)
       endif
    end if

! -----------------------------------------------------------------
! Calcul niveaux pression modele
! -----------------------------------------------------------------

!    ....  Calcul de pls , pression au milieu des couches ,en Pascals
    IF (guide_plevs.EQ.1) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
            DO j=jjbu,jjeu
                DO i =1, iip1
                    pls(i,j,l)=(ap(l)+ap(l+1))/2.+psi(i,j)*(bp(l)+bp(l+1))/2.
                ENDDO
            ENDDO
        ENDDO
    ELSE
        unskap=1./kappa
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
   DO l = 1, llm
       DO j=jjbu,jjeu
           DO i =1, iip1
               pls(i,j,l) = preff * ( pk(i,j,l)/cpp) ** unskap
           ENDDO
       ENDDO
   ENDDO
    ENDIF

!   calcul des pressions pour les grilles u et v
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    do l=1,llm
        do j=jjbu,jjeu
            do i=1,iip1
                pext(i,j,l)=pls(i,j,l)*aire(i,j)
            enddo
        enddo
    enddo

     CALL Register_Hallo_u(pext,llm,1,2,2,1,Req)
     CALL SendRequest(Req)
!$OMP BARRIER
     CALL WaitRequest(Req)
!$OMP BARRIER

    call massbar_loc(pext, pbarx, pbary )
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    do l=1,llm
        do j=jjbu,jjeu
            do i=1,iip1
                plunc(i,j,l)=pbarx(i,j,l)/aireu(i,j)
                plsnc(i,j,l)=pls(i,j,l)
            enddo
        enddo
    enddo
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
    do l=1,llm
        do j=jjbv,jjev
            do i=1,iip1
                plvnc(i,j,l)=pbary(i,j,l)/airev(i,j)
            enddo
        enddo
    enddo

! -----------------------------------------------------------------
! Interpolation verticale champs guidage sur niveaux modele
! Conversion en variables gcm (ucov, vcov...)
! -----------------------------------------------------------------
    if (guide_P) then
!$OMP MASTER
        do j=jjbu,jjeu
            do i=1,iim
                ij=(j-1)*iip1+i
                psgui1(ij)=psnat1(i,j)
                psgui2(ij)=psnat2(i,j)
            enddo
            psgui1(iip1*j)=psnat1(1,j)
            psgui2(iip1*j)=psnat2(1,j)
        enddo
!$OMP END MASTER
!$OMP BARRIER
    endif

    IF (guide_T) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*psnat2(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*psnat1(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iip1
                        plnc2(i,j,l)=pnat2(i,j,l)
                        plnc1(i,j,l)=pnat1(i,j,l)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        ! Interpolation verticale
!$OMP MASTER
        CALL pres2lev(tnat1(:,jjbu:jjeu,:),zu1(:,jjbu:jjeu,:),nlevnc,llm,           &
                    plnc1(:,jjbu:jjeu,:),plsnc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
        CALL pres2lev(tnat2(:,jjbu:jjeu,:),zu2(:,jjbu:jjeu,:),nlevnc,llm,           &
                    plnc2(:,jjbu:jjeu,:),plsnc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
!$OMP END MASTER
!$OMP BARRIER
        ! Conversion en variables GCM
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
            do j=jjbu,jjeu
                IF (guide_teta) THEN
                    do i=1,iim
                        ij=(j-1)*iip1+i
                        tgui1(ij,l)=zu1(i,j,l)
                        tgui2(ij,l)=zu2(i,j,l)
                    enddo
                ELSE
                    do i=1,iim
                        ij=(j-1)*iip1+i
                        tgui1(ij,l)=zu1(i,j,l)*cpp/pk(i,j,l)
                        tgui2(ij,l)=zu2(i,j,l)*cpp/pk(i,j,l)
                    enddo
                ENDIF
                tgui1(j*iip1,l)=tgui1((j-1)*iip1+1,l)
                tgui2(j*iip1,l)=tgui2((j-1)*iip1+1,l)
            enddo
            if (pole_nord) then
              do i=1,iip1
                tgui1(i,l)=tgui1(1,l)
                tgui2(i,l)=tgui2(1,l)
              enddo
            endif
            if (pole_sud) then
              do i=1,iip1
                tgui1(ip1jm+i,l)=tgui1(ip1jm+1,l)
                tgui2(ip1jm+i,l)=tgui2(ip1jm+1,l)
              enddo
           endif
        enddo
    ENDIF

    IF (guide_Q) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*psnat2(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*psnat1(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iip1
                        plnc2(i,j,l)=pnat2(i,j,l)
                        plnc1(i,j,l)=pnat1(i,j,l)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        ! Interpolation verticale
!$OMP MASTER
        CALL pres2lev(qnat1(:,jjbu:jjeu,:),zu1(:,jjbu:jjeu,:),nlevnc,llm,             &
                      plnc1(:,jjbu:jjeu,:),plsnc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
        CALL pres2lev(qnat2(:,jjbu:jjeu,:),zu2(:,jjbu:jjeu,:),nlevnc,llm,             &
                      plnc2(:,jjbu:jjeu,:),plsnc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
!$OMP END MASTER
!$OMP BARRIER

        ! Conversion en variables GCM
        ! On suppose qu'on a la bonne variable dans le fichier de guidage:
        ! Hum.Rel si guide_hr, Hum.Spec. sinon.
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
            do j=jjbu,jjeu
                do i=1,iim
                    ij=(j-1)*iip1+i
                    qgui1(ij,l)=zu1(i,j,l)
                    qgui2(ij,l)=zu2(i,j,l)
                enddo
                qgui1(j*iip1,l)=qgui1((j-1)*iip1+1,l)
                qgui2(j*iip1,l)=qgui2((j-1)*iip1+1,l)
            enddo
            if (pole_nord) then
              do i=1,iip1
                qgui1(i,l)=qgui1(1,l)
                qgui2(i,l)=qgui2(1,l)
              enddo
            endif
            if (pole_sud) then
              do i=1,iip1
                qgui1(ip1jm+i,l)=qgui1(ip1jm+1,l)
                qgui2(ip1jm+i,l)=qgui2(ip1jm+1,l)
              enddo
            endif
        enddo
        IF (guide_hr) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          do l=1,llm
            CALL q_sat(iip1*jjnu,teta(:,jjbu:jjeu,l)*pk(:,jjbu:jjeu,l)/cpp,       &
                       plsnc(:,jjbu:jjeu,l),qsat(ijbu:ijeu,l))
            qgui1(ijbu:ijeu,l)=qgui1(ijbu:ijeu,l)*qsat(ijbu:ijeu,l)*0.01 !hum. rel. en %
            qgui2(ijbu:ijeu,l)=qgui2(ijbu:ijeu,l)*qsat(ijbu:ijeu,l)*0.01
          enddo

        ENDIF
    ENDIF

    IF (guide_u) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iim
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*(psnat2(i,j)*aire(i,j)*alpha1p2(i,j) &
                       &           +psnat2(i+1,j)*aire(i+1,j)*alpha3p4(i+1,j))/aireu(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*(psnat1(i,j)*aire(i,j)*alpha1p2(i,j) &
                       &           +psnat1(i+1,j)*aire(i+1,j)*alpha3p4(i+1,j))/aireu(i,j)
                    ENDDO
                    plnc2(iip1,j,l)=plnc2(1,j,l)
                    plnc1(iip1,j,l)=plnc1(1,j,l)
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbu,jjeu
                    DO i=1,iim
                        plnc2(i,j,l)=(pnat2(i,j,l)*aire(i,j)*alpha1p2(i,j) &
                       & +pnat2(i+1,j,l)*aire(i,j)*alpha3p4(i+1,j))/aireu(i,j)
                        plnc1(i,j,l)=(pnat1(i,j,l)*aire(i,j)*alpha1p2(i,j) &
                       & +pnat1(i+1,j,l)*aire(i,j)*alpha3p4(i+1,j))/aireu(i,j)
                    ENDDO
                    plnc2(iip1,j,l)=plnc2(1,j,l)
                    plnc1(iip1,j,l)=plnc1(1,j,l)
                ENDDO
            ENDDO
        ENDIF

        ! Interpolation verticale
!$OMP MASTER
        CALL pres2lev(unat1(:,jjbu:jjeu,:),zu1(:,jjbu:jjeu,:),nlevnc,llm,            &
                      plnc1(:,jjbu:jjeu,:),plunc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
        CALL pres2lev(unat2(:,jjbu:jjeu,:),zu2(:,jjbu:jjeu,:),nlevnc,llm,            &
                      plnc2(:,jjbu:jjeu,:),plunc(:,jjbu:jjeu,:),iip1,jjnu,invert_p)
!$OMP END MASTER
!$OMP BARRIER

        ! Conversion en variables GCM
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
            do j=jjbu,jjeu
                do i=1,iim
                    ij=(j-1)*iip1+i
                    ugui1(ij,l)=zu1(i,j,l)*cu(i,j)
                    ugui2(ij,l)=zu2(i,j,l)*cu(i,j)
                enddo
                ugui1(j*iip1,l)=ugui1((j-1)*iip1+1,l)
                ugui2(j*iip1,l)=ugui2((j-1)*iip1+1,l)
            enddo
            if (pole_nord) then
              do i=1,iip1
                ugui1(i,l)=0.
                ugui2(i,l)=0.
              enddo
            endif
            if (pole_sud) then
              do i=1,iip1
                ugui1(ip1jm+i,l)=0.
                ugui2(ip1jm+i,l)=0.
              enddo
            endif
        enddo
    ENDIF

    IF (guide_v) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
         CALL Register_Hallo_u(psnat1,1,1,2,2,1,Req)
         CALL Register_Hallo_u(psnat2,1,1,2,2,1,Req)
         CALL SendRequest(Req)
!$OMP BARRIER
         CALL WaitRequest(Req)
!$OMP BARRIER
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbv,jjev
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*(psnat2(i,j)*aire(i,j)*alpha2p3(i,j) &
                       &           +psnat2(i,j+1)*aire(i,j+1)*alpha1p4(i,j+1))/airev(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*(psnat1(i,j)*aire(i,j)*alpha2p3(i,j) &
                       &           +psnat1(i,j+1)*aire(i,j+1)*alpha1p4(i,j+1))/airev(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
         CALL Register_Hallo_u(pnat1,llm,1,2,2,1,Req)
         CALL Register_Hallo_u(pnat2,llm,1,2,2,1,Req)
         CALL SendRequest(Req)
!$OMP BARRIER
         CALL WaitRequest(Req)
!$OMP BARRIER
!$OMP DO
            DO l=1,nlevnc
                DO j=jjbv,jjev
                    DO i=1,iip1
                        plnc2(i,j,l)=(pnat2(i,j,l)*aire(i,j)*alpha2p3(i,j) &
                       & +pnat2(i,j+1,l)*aire(i,j)*alpha1p4(i,j+1))/airev(i,j)
                        plnc1(i,j,l)=(pnat1(i,j,l)*aire(i,j)*alpha2p3(i,j) &
                       & +pnat1(i,j+1,l)*aire(i,j)*alpha1p4(i,j+1))/airev(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
        ! Interpolation verticale

!$OMP MASTER
        CALL pres2lev(vnat1(:,jjbv:jjev,:),zv1(:,jjbv:jjev,:),nlevnc,llm,             &
                      plnc1(:,jjbv:jjev,:),plvnc(:,jjbv:jjev,:),iip1,jjnv,invert_p)
        CALL pres2lev(vnat2(:,jjbv:jjev,:),zv2(:,jjbv:jjev,:),nlevnc,llm,             &
                      plnc2(:,jjbv:jjev,:),plvnc(:,jjbv:jjev,:),iip1,jjnv,invert_p)
!$OMP END MASTER
!$OMP BARRIER
        ! Conversion en variables GCM
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
            do j=jjbv,jjev
                do i=1,iim
                    ij=(j-1)*iip1+i
                    vgui1(ij,l)=zv1(i,j,l)*cv(i,j)
                    vgui2(ij,l)=zv2(i,j,l)*cv(i,j)
                enddo
                vgui1(j*iip1,l)=vgui1((j-1)*iip1+1,l)
                vgui2(j*iip1,l)=vgui2((j-1)*iip1+1,l)
            enddo
        enddo
    ENDIF


  END SUBROUTINE guide_interp

!=======================================================================
  SUBROUTINE tau2alpha(typ,pim,jjb,jje,factt,taumin,taumax,alpha)

! Calcul des constantes de rappel alpha (=1/tau)

    use comconst_mod, only: pi
    use serre_mod, only: clat, clon, grossismx, grossismy

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
    USE comgeom2_mod_h
    implicit none


! input arguments :
    INTEGER, INTENT(IN) :: typ    ! u(2),v(3), ou scalaire(1)
    INTEGER, INTENT(IN) :: pim ! dimensions en lon
    INTEGER, INTENT(IN) :: jjb,jje ! dimensions en lat
    REAL, INTENT(IN)    :: factt   ! pas de temps en fraction de jour
    REAL, INTENT(IN)    :: taumin,taumax
! output arguments:
    REAL, DIMENSION(pim,jjb:jje), INTENT(OUT) :: alpha

!  local variables:
    LOGICAL, SAVE               :: first=.TRUE.
    REAL, SAVE                  :: gamma,dxdy_min,dxdy_max
    REAL, DIMENSION (iip1,jjp1) :: zdx,zdy
    REAL, DIMENSION (iip1,jjp1) :: dxdys,dxdyu
    REAL, DIMENSION (iip1,jjm)  :: dxdyv
    real dxdy_
    real zlat,zlon
    real alphamin,alphamax,xi
    integer i,j,ilon,ilat
    character(len=20),parameter :: modname="tau2alpha"


    alphamin=factt/taumax
    alphamax=factt/taumin
    IF (guide_reg.OR.guide_add) THEN
        alpha=alphamax
!-----------------------------------------------------------------------
! guide_reg: alpha=alpha_min dans region, 0. sinon.
!-----------------------------------------------------------------------
        IF (guide_reg) THEN
            do j=jjb,jje
                do i=1,pim
                    if (typ.eq.2) then
                       zlat=rlatu(j)*180./pi
                       zlon=rlonu(i)*180./pi
                    elseif (typ.eq.1) then
                       zlat=rlatu(j)*180./pi
                       zlon=rlonv(i)*180./pi
                    elseif (typ.eq.3) then
                       zlat=rlatv(j)*180./pi
                       zlon=rlonv(i)*180./pi
                    endif
                    alpha(i,j)=alphamax/16.* &
                              (1.+tanh((zlat-lat_min_g)/tau_lat))* &
                              (1.+tanh((lat_max_g-zlat)/tau_lat))* &
                              (1.+tanh((zlon-lon_min_g)/tau_lon))* &
                              (1.+tanh((lon_max_g-zlon)/tau_lon))
                enddo
            enddo
        ENDIF
    ELSE
!-----------------------------------------------------------------------
! Sinon, alpha varie entre alpha_min et alpha_max suivant le zoom.
!-----------------------------------------------------------------------
!Calcul de l'aire des mailles
        do j=2,jjm
            do i=2,iip1
               zdx(i,j)=0.5*(cu(i-1,j)+cu(i,j))/cos(rlatu(j))
            enddo
            zdx(1,j)=zdx(iip1,j)
        enddo
        do j=2,jjm
            do i=1,iip1
               zdy(i,j)=0.5*(cv(i,j-1)+cv(i,j))
            enddo
        enddo
        do i=1,iip1
            zdx(i,1)=zdx(i,2)
            zdx(i,jjp1)=zdx(i,jjm)
            zdy(i,1)=zdy(i,2)
            zdy(i,jjp1)=zdy(i,jjm)
        enddo
        do j=1,jjp1
            do i=1,iip1
               dxdys(i,j)=sqrt(zdx(i,j)*zdx(i,j)+zdy(i,j)*zdy(i,j))
            enddo
        enddo
        IF (typ.EQ.2) THEN
            do j=1,jjp1
                do i=1,iim
                   dxdyu(i,j)=0.5*(dxdys(i,j)+dxdys(i+1,j))
                enddo
                dxdyu(iip1,j)=dxdyu(1,j)
            enddo
        ENDIF
        IF (typ.EQ.3) THEN
            do j=1,jjm
                do i=1,iip1
                   dxdyv(i,j)=0.5*(dxdys(i,j)+dxdys(i,j+1))
                enddo
            enddo
        ENDIF
! Premier appel: calcul des aires min et max et de gamma.
        IF (first) THEN
            first=.FALSE.
            ! coordonnees du centre du zoom
            CALL coordij(clon,clat,ilon,ilat)
            ! aire de la maille au centre du zoom
            dxdy_min=dxdys(ilon,ilat)
            ! dxdy maximale de la maille
            dxdy_max=0.
            do j=1,jjp1
                do i=1,iip1
                     dxdy_max=max(dxdy_max,dxdys(i,j))
                enddo
            enddo
            ! Calcul de gamma
            if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
              write(*,*)trim(modname)//' ATTENTION modele peu zoome'
              write(*,*)trim(modname)//' ATTENTION on prend une constante de guidage cste'
              gamma=0.
            else
              gamma=(dxdy_max-2.*dxdy_min)/(dxdy_max-dxdy_min)
              write(*,*)trim(modname)//' gamma=',gamma
              if (gamma.lt.1.e-5) then
                write(*,*)trim(modname)//' gamma =',gamma,'<1e-5'
                CALL abort_gcm("guide_loc_mod","stopped",1)
              endif
              gamma=log(0.5)/log(gamma)
              if (gamma4) then
                gamma=min(gamma,4.)
              endif
              write(*,*)trim(modname)//' gamma=',gamma
            endif
        ENDIF !first

        do j=jjb,jje
            do i=1,pim
                if (typ.eq.1) then
                   dxdy_=dxdys(i,j)
                   zlat=rlatu(j)*180./pi
                elseif (typ.eq.2) then
                   dxdy_=dxdyu(i,j)
                   zlat=rlatu(j)*180./pi
                elseif (typ.eq.3) then
                   dxdy_=dxdyv(i,j)
                   zlat=rlatv(j)*180./pi
                endif
                if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
                ! pour une grille reguliere, xi=xxx**0=1 -> alpha=alphamin
                    alpha(i,j)=alphamin
                else
                    xi=((dxdy_max-dxdy_)/(dxdy_max-dxdy_min))**gamma
                    xi=min(xi,1.)
                    if(lat_min_g.le.zlat .and. zlat.le.lat_max_g) then
                        alpha(i,j)=xi*alphamin+(1.-xi)*alphamax
                    else
                        alpha(i,j)=0.
                    endif
                endif
            enddo
        enddo
    ENDIF ! guide_reg

    if (.not. guide_add) alpha = 1. - exp(- alpha)

  END SUBROUTINE tau2alpha

!=======================================================================
  SUBROUTINE guide_read(timestep)
    USE netcdf, ONLY: nf90_get_var
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
    USE pyfld, ONLY: py_du_caldyn, py_dv_caldyn
    USE pycpl, ONLY: lk_pycpl
IMPLICIT NONE



    INTEGER, INTENT(IN)   :: timestep

    LOGICAL, SAVE         :: first=.TRUE.
! Identification fichiers et variables NetCDF:
    INTEGER, SAVE         :: ncidu,varidu,ncidv,varidv,ncidp,varidp
    INTEGER, SAVE         :: ncidQ,varidQ,ncidt,varidt,ncidps,varidps
    INTEGER               :: ncidpl,varidpl,varidap,varidbp,dimid,lendim
! Variables auxiliaires NetCDF:
    INTEGER, DIMENSION(4) :: start,count
    INTEGER               :: status,rcode
    INTEGER               :: j,l
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_read'
    CHARACTER (len = 20)   :: namedim
    abort_message='pb in guide_read'

! -----------------------------------------------------------------
! Premier appel: initialisation de la lecture des fichiers
! -----------------------------------------------------------------
    if (first) then
         ncidpl=-99
         write(*,*) trim(modname)//': opening nudging files '
! Ap et Bp si Niveaux de pression hybrides
         if (guide_plevs.EQ.1) then
             write(*,*) trim(modname)//' Reading nudging on model levels'
             rcode = nf90_open('apbp.nc', nf90_nowrite, ncidpl)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'AP', varidap)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no AP variable in file apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'BP', varidbp)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no BP variable in file apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidpl,varidap',ncidpl,varidap
         endif

! Pression si guidage sur niveaux P variables
         if (guide_plevs.EQ.2) then
             rcode = nf90_open('P.nc', nf90_nowrite, ncidp)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidp, 'PRES', varidp)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no PRES variable in file P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidp,varidp',ncidp,varidp
             if (ncidpl.eq.-99) ncidpl=ncidp
         endif

! Vent zonal
         if (guide_u) then
             rcode = nf90_open('u.nc', nf90_nowrite, ncidu)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidu, 'UWND', varidu)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no UWND variable in file u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidu,varidu',ncidu,varidu
             if (ncidpl.eq.-99) ncidpl=ncidu


             status=NF90_INQ_DIMID(ncidu, "LONU", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidu,dimid,namedim,lendim)
             IF (lendim .NE. iip1) THEN
                abort_message='dimension LONU different from iip1 in u.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

             status=NF90_INQ_DIMID(ncidu, "LATU", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidu,dimid,namedim,lendim)
             IF (lendim .NE. jjp1) THEN
                abort_message='dimension LATU different from jjp1 in u.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

         endif

! Vent meridien
         if (guide_v) then
             rcode = nf90_open('v.nc', nf90_nowrite, ncidv)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidv, 'VWND', varidv)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no VWND variable in file v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidv,varidv',ncidv,varidv
             if (ncidpl.eq.-99) ncidpl=ncidv

             status=NF90_INQ_DIMID(ncidv, "LONV", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidv,dimid,namedim,lendim)

                IF (lendim .NE. iip1) THEN
                abort_message='dimension LONV different from iip1 in v.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF


             status=NF90_INQ_DIMID(ncidv, "LATV", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidv,dimid,namedim,lendim)
             IF (lendim .NE. jjm) THEN
                abort_message='dimension LATV different from jjm in v.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

        endif

! Temperature
         if (guide_T) then
             rcode = nf90_open('T.nc', nf90_nowrite, ncidt)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidt, 'AIR', varidt)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no AIR variable in file T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidT,varidT',ncidt,varidt
             if (ncidpl.eq.-99) ncidpl=ncidt

             status=NF90_INQ_DIMID(ncidt, "LONV", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidt,dimid,namedim,lendim)
             IF (lendim .NE. iip1) THEN
                abort_message='dimension LONV different from iip1 in T.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

             status=NF90_INQ_DIMID(ncidt, "LATU", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidt,dimid,namedim,lendim)
             IF (lendim .NE. jjp1) THEN
                abort_message='dimension LATU different from jjp1 in T.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

         endif

! Humidite
         if (guide_Q) then
             rcode = nf90_open('hur.nc', nf90_nowrite, ncidQ)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidQ, 'RH', varidQ)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no RH variable in file hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidQ,varidQ',ncidQ,varidQ
             if (ncidpl.eq.-99) ncidpl=ncidQ


             status=NF90_INQ_DIMID(ncidQ, "LONV", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidQ,dimid,namedim,lendim)
             IF (lendim .NE. iip1) THEN
                abort_message='dimension LONV different from iip1 in hur.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF

             status=NF90_INQ_DIMID(ncidQ, "LATU", dimid)
             status=NF90_INQUIRE_DIMENSION(ncidQ,dimid,namedim,lendim)
             IF (lendim .NE. jjp1) THEN
                abort_message='dimension LATU different from jjp1 in hur.nc'
                CALL abort_gcm(modname,abort_message,1)
             ENDIF


         endif
! Pression de surface
         if ((guide_P).OR.(guide_plevs.EQ.1)) then
             rcode = nf90_open('ps.nc', nf90_nowrite, ncidps)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no file ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidps, 'SP', varidps)
             IF (rcode.NE.nf90_noerr) THEN
              abort_message='Nudging: error -> no SP variable in file ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             write(*,*) trim(modname)//' ncidps,varidps',ncidps,varidps
         endif
! Coordonnee verticale
         if (guide_plevs.EQ.0) then
              rcode = nf90_inq_varid(ncidpl, 'LEVEL', varidpl)
              IF (rcode.NE.0) rcode = nf90_inq_varid(ncidpl, 'PRESSURE', varidpl)
              write(*,*) trim(modname)//' ncidpl,varidpl',ncidpl,varidpl
         endif
! Coefs ap, bp pour calcul de la pression aux differents niveaux
!function nf90_get_var(ncid, varid, values, start, count, stride, map)
!  integer,                         intent( in) :: ncid, varid
!  any valid type, scalar or array of any rank, &
!                                   intent(out) :: values
!  integer, dimension(:), optional, intent( in) :: start, count, stride, map
!  integer                                      :: nf90_get_var
         IF (guide_plevs.EQ.1) THEN
             status = nf90_get_var(ncidpl, varidap, apnc, [1], [nlevnc])
             status = nf90_get_var(ncidpl, varidbp, bpnc, [1], [nlevnc])
         ELSEIF (guide_plevs.EQ.0) THEN
             status = nf90_get_var(ncidpl, varidpl, apnc, [1], [nlevnc])
!FC Pour les corrections la pression est deja en Pascals on commente la ligne ci-dessous
             IF(convert_Pa) apnc=apnc*100.! conversion en Pascals
             bpnc(:)=0.
         ENDIF
         first=.FALSE.
     ENDIF ! (first)

! -----------------------------------------------------------------
!   lecture des champs u, v, T, Q, ps
! -----------------------------------------------------------------

!  dimensions pour les champs scalaires et le vent zonal
     start(1)=1
     start(2)=jjb_u
     start(3)=1
     start(4)=timestep

     count(1)=iip1
     count(2)=jjnb_u
     count(3)=nlevnc
     count(4)=1

     IF (invert_y) start(2)=jjp1-jje_u+1
! Pression
     if (guide_plevs.EQ.2) then
         status = nf90_get_var(ncidp, varidp, pnat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,pnat2)
         ENDIF
     endif

!  Vent zonal
     if (guide_u) then
         status = nf90_get_var(ncidu, varidu, unat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
            CALL invert_lat(iip1,jjnb_u,nlevnc,unat2)
         ENDIF
         IF (lk_pycpl .AND. ALLOCATED(py_du_caldyn)) THEN
             ! +++ Python coupling: nudging target = python field
             DO l = 1, llm
                DO j = jjb_u, jje_u
                   unat2(:,j,l) = py_du_caldyn((j-1)*iip1+1 : j*iip1, l)
                ENDDO
             ENDDO
         ENDIF
     endif


!  Temperature
     if (guide_T) then
         status = nf90_get_var(ncidt, varidt, tnat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,tnat2)
         ENDIF
     endif

!  Humidite
     if (guide_Q) then
         status = nf90_get_var(ncidQ, varidQ, qnat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,qnat2)
         ENDIF

     endif

!  Vent meridien
     if (guide_v) then
         start(2)=jjb_v
         count(2)=jjnb_v
         IF (invert_y) start(2)=jjm-jje_v+1
         status = nf90_get_var(ncidv, varidv, vnat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
            CALL invert_lat(iip1,jjnb_v,nlevnc,vnat2)
         ENDIF
         IF (lk_pycpl .AND. ALLOCATED(py_dv_caldyn)) THEN
            ! +++ Python coupling: nudging target = python field
             DO l = 1, llm
                DO j = jjb_v, jje_v
                   vnat2(:,j,l) = py_dv_caldyn((j-1)*iip1+1 : j*iip1, l)
                ENDDO
             ENDDO
         ENDIF
     endif

!  Pression de surface
     if ((guide_P).OR.(guide_plevs.EQ.1))  then
         start(2)=jjb_u
         start(3)=timestep
         start(4)=0
         count(2)=jjnb_u
         count(3)=1
         count(4)=0
         IF (invert_y) start(2)=jjp1-jje_u+1
         status = nf90_get_var(ncidps, varidps, psnat2, start, count)
         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,1,psnat2)
         ENDIF
     endif

  END SUBROUTINE guide_read

!=======================================================================
  SUBROUTINE guide_read2D(timestep)
    USE netcdf, ONLY: nf90_get_var
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    INTEGER, INTENT(IN)   :: timestep

    LOGICAL, SAVE         :: first=.TRUE.
! Identification fichiers et variables NetCDF:
    INTEGER, SAVE         :: ncidu,varidu,ncidv,varidv,ncidp,varidp
    INTEGER, SAVE         :: ncidQ,varidQ,ncidt,varidt,ncidps,varidps
    INTEGER               :: ncidpl,varidpl,varidap,varidbp
! Variables auxiliaires NetCDF:
    INTEGER, DIMENSION(4) :: start,count
    INTEGER               :: status,rcode
! Variables for 3D extension:
    REAL, DIMENSION (jjb_u:jje_u,llm)  :: zu
    REAL, DIMENSION (jjb_v:jje_v,llm)  :: zv
    INTEGER               :: i
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_read2D'
    abort_message='pb in guide_read2D'

! -----------------------------------------------------------------
! Premier appel: initialisation de la lecture des fichiers
! -----------------------------------------------------------------
    if (first) then
         ncidpl=-99
         write(*,*)trim(modname)//' : opening nudging files '
! Ap et Bp si niveaux de pression hybrides
         if (guide_plevs.EQ.1) then
           write(*,*)trim(modname)//' Reading nudging on model levels'
           rcode = nf90_open('apbp.nc', nf90_nowrite, ncidpl)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file apbp.nc'
           CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidpl, 'AP', varidap)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no AP variable in file apbp.nc'
           CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidpl, 'BP', varidbp)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no BP variable in file apbp.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//'ncidpl,varidap',ncidpl,varidap
         endif
! Pression
         if (guide_plevs.EQ.2) then
           rcode = nf90_open('P.nc', nf90_nowrite, ncidp)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file P.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidp, 'PRES', varidp)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no PRES variable in file P.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidp,varidp',ncidp,varidp
           if (ncidpl.eq.-99) ncidpl=ncidp
         endif
! Vent zonal
         if (guide_u) then
           rcode = nf90_open('u.nc', nf90_nowrite, ncidu)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file u.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidu, 'UWND', varidu)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no UWND variable in file u.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidu,varidu',ncidu,varidu
           if (ncidpl.eq.-99) ncidpl=ncidu
         endif

! Vent meridien
         if (guide_v) then
           rcode = nf90_open('v.nc', nf90_nowrite, ncidv)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file v.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidv, 'VWND', varidv)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no VWND variable in file v.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidv,varidv',ncidv,varidv
           if (ncidpl.eq.-99) ncidpl=ncidv
        endif
! Temperature
         if (guide_T) then
           rcode = nf90_open('T.nc', nf90_nowrite, ncidt)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file T.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidt, 'AIR', varidt)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no AIR variable in file T.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidT,varidT',ncidt,varidt
           if (ncidpl.eq.-99) ncidpl=ncidt
         endif
! Humidite
         if (guide_Q) then
           rcode = nf90_open('hur.nc', nf90_nowrite, ncidQ)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file hur.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidQ, 'RH', varidQ)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no RH,variable in file hur.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidQ,varidQ',ncidQ,varidQ
           if (ncidpl.eq.-99) ncidpl=ncidQ
         endif
! Pression de surface
         if ((guide_P).OR.(guide_plevs.EQ.1)) then
           rcode = nf90_open('ps.nc', nf90_nowrite, ncidps)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no file ps.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           rcode = nf90_inq_varid(ncidps, 'SP', varidps)
           IF (rcode.NE.nf90_noerr) THEN
             abort_message='Nudging: error -> no SP variable in file ps.nc'
             CALL abort_gcm(modname,abort_message,1)
           ENDIF
           write(*,*)trim(modname)//' ncidps,varidps',ncidps,varidps
         endif
! Coordonnee verticale
         if (guide_plevs.EQ.0) then
           rcode = nf90_inq_varid(ncidpl, 'LEVEL', varidpl)
           IF (rcode.NE.0) rcode = nf90_inq_varid(ncidpl, 'PRESSURE', varidpl)
           write(*,*)trim(modname)//' ncidpl,varidpl',ncidpl,varidpl
         endif
! Coefs ap, bp pour calcul de la pression aux differents niveaux
         if (guide_plevs.EQ.1) then
             status = nf90_get_var(ncidpl, varidap, apnc, [1], [nlevnc])
             status = nf90_get_var(ncidpl, varidbp, bpnc, [1], [nlevnc])
         elseif (guide_plevs.EQ.0) THEN
             status = nf90_get_var(ncidpl, varidpl, apnc, [1], [nlevnc])
             apnc=apnc*100.! conversion en Pascals
             bpnc(:)=0.
         endif
         first=.FALSE.
     endif ! (first)

! -----------------------------------------------------------------
!   lecture des champs u, v, T, Q, ps
! -----------------------------------------------------------------

!  dimensions pour les champs scalaires et le vent zonal
     start(1)=1
     start(2)=jjb_u
     start(3)=1
     start(4)=timestep

     count(1)=1
     count(2)=jjnb_u
     count(3)=nlevnc
     count(4)=1

     IF (invert_y) start(2)=jjp1-jje_u+1
!  Pression
     if (guide_plevs.EQ.2) then
         status = nf90_get_var(ncidp, varidp, zu, start, count)
         DO i=1,iip1
             pnat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,pnat2)
         ENDIF
     endif
!  Vent zonal
     if (guide_u) then
         status = nf90_get_var(ncidu, varidu, zu, start, count)
         DO i=1,iip1
             unat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,unat2)
         ENDIF
     endif


!  Temperature
     if (guide_T) then
         status = nf90_get_var(ncidt, varidt, zu, start, count)
         DO i=1,iip1
             tnat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,tnat2)
         ENDIF
     endif

!  Humidite
     if (guide_Q) then
         status = nf90_get_var(ncidQ, varidQ, zu, start, count)
         DO i=1,iip1
             qnat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,nlevnc,qnat2)
         ENDIF
     endif

!  Vent meridien
     if (guide_v) then
         start(2)=jjb_v
         count(2)=jjnb_v
         IF (invert_y) start(2)=jjm-jje_v+1
         status = nf90_get_var(ncidv, varidv, zv, start, count)
         DO i=1,iip1
             vnat2(i,:,:)=zv(:,:)
         ENDDO

         IF (invert_y) THEN

!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_v,nlevnc,vnat2)
         ENDIF
     endif

!  Pression de surface
     if ((guide_P).OR.(guide_plevs.EQ.1))  then
         start(2)=jjb_u
         start(3)=timestep
         start(4)=0
         count(2)=jjnb_u
         count(3)=1
         count(4)=0
         IF (invert_y) start(2)=jjp1-jje_u+1
         status = nf90_get_var(ncidps, varidps, zu(:, 1), start, count)
         DO i=1,iip1
             psnat2(i,:)=zu(:,1)
         ENDDO

         IF (invert_y) THEN
!           PRINT*,"Invertion impossible actuellement"
!           CALL abort_gcm(modname,abort_message,1)
           CALL invert_lat(iip1,jjnb_u,1,psnat2)
         ENDIF
     endif

  END SUBROUTINE guide_read2D

!=======================================================================
  SUBROUTINE guide_out(varname,hsize,vsize,field_loc,factt)
    USE parallel_lmdz
    USE mod_hallo, ONLY : gather_field_u, gather_field_v
    USE comconst_mod, ONLY: pi
    USE comvert_mod, ONLY: presnivs
    use netcdf95, only: nf95_def_var, nf95_put_var
    use netcdf, only: nf90_float, nf90_put_var

    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
    USE comgeom2_mod_h
IMPLICIT NONE




    ! Variables entree
    CHARACTER*(*), INTENT(IN)                      :: varname
    INTEGER,   INTENT (IN)                         :: hsize,vsize
!   REAL, DIMENSION (iip1,hsize,vsize), INTENT(IN) :: field_loc
    REAL, DIMENSION (:,:), INTENT(IN) :: field_loc
    REAL factt

    ! Variables locales
    INTEGER, SAVE :: timestep=0
    ! Identites fichier netcdf
    INTEGER       :: nid, id_lonu, id_lonv, id_latu, id_latv, id_tim, id_lev
    INTEGER       :: vid_lonu,vid_lonv,vid_latu,vid_latv,vid_cu,vid_cv,vid_lev
    INTEGER       :: vid_au,vid_av, varid_alpha_t, varid_alpha_q
    INTEGER, DIMENSION (3) :: dim3
    INTEGER, DIMENSION (4) :: dim4,count,start
    INTEGER                :: ierr, varid,l
    REAL zu(ip1jmp1),zv(ip1jm), zt(iip1, jjp1), zq(iip1, jjp1)
    REAL, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: field_glo
    CHARACTER(LEN=20),PARAMETER :: modname="guide_out"

!$OMP MASTER
    ALLOCATE(field_glo(iip1,hsize,vsize))
!$OMP END MASTER
!$OMP BARRIER

!    write(*,*)trim(modname)//' after allocation ',hsize,vsize

    IF (hsize==jjp1) THEN
        CALL gather_field_u(field_loc,field_glo,vsize)
    ELSE IF (hsize==jjm) THEN
       CALL gather_field_v(field_loc,field_glo, vsize)
    ENDIF

!    write(*,*)trim(modname)//' after gather '
    CALL Gather_field_u(alpha_u,zu,1)
    CALL Gather_field_u(alpha_t,zt,1)
    CALL Gather_field_u(alpha_q,zq,1)
    CALL Gather_field_v(alpha_v,zv,1)

    IF (mpi_rank >  0) THEN
!$OMP MASTER
       DEALLOCATE(field_glo)
!$OMP END MASTER
!$OMP BARRIER

       RETURN
    ENDIF

!$OMP MASTER
    IF (timestep.EQ.0) THEN
! ----------------------------------------------
! initialisation fichier de sortie
! ----------------------------------------------
! Ouverture du fichier
        ierr=nf90_create("guide_ins.nc",IOR(nf90_clobber,nf90_64bit_offset),nid)
! Definition des dimensions
        ierr=nf90_def_dim(nid,"LONU",iip1,id_lonu)
        ierr=nf90_def_dim(nid,"LONV",iip1,id_lonv)
        ierr=nf90_def_dim(nid,"LATU",jjp1,id_latu)
        ierr=nf90_def_dim(nid,"LATV",jjm,id_latv)
        ierr=nf90_def_dim(nid,"LEVEL",llm,id_lev)
        ierr=nf90_def_dim(nid,"TIME",nf90_unlimited,id_tim)

! Creation des variables dimensions
        ierr=nf90_def_var(nid,"LONU",nf90_float,id_lonu,vid_lonu)
        ierr=nf90_def_var(nid,"LONV",nf90_float,id_lonv,vid_lonv)
        ierr=nf90_def_var(nid,"LATU",nf90_float,id_latu,vid_latu)
        ierr=nf90_def_var(nid,"LATV",nf90_float,id_latv,vid_latv)
        ierr=nf90_def_var(nid,"LEVEL",nf90_float,id_lev,vid_lev)
        ierr=nf90_def_var(nid,"cu",nf90_float,(/id_lonu,id_latu/),vid_cu)
        ierr=nf90_def_var(nid,"cv",nf90_float,(/id_lonv,id_latv/),vid_cv)
        ierr=nf90_def_var(nid,"au",nf90_float,(/id_lonu,id_latu/),vid_au)
        ierr=nf90_def_var(nid,"av",nf90_float,(/id_lonv,id_latv/),vid_av)
        call nf95_def_var(nid, "alpha_T", nf90_float, (/id_lonv, id_latu/), &
             varid_alpha_t)
        call nf95_def_var(nid, "alpha_q", nf90_float, (/id_lonv, id_latu/), &
             varid_alpha_q)

        ierr=nf90_enddef(nid)

! Enregistrement des variables dimensions
        ierr = nf90_put_var(nid, vid_lonu, rlonu * 180. / pi)
        ierr = nf90_put_var(nid, vid_lonv, rlonv * 180. / pi)
        ierr = nf90_put_var(nid, vid_latu, rlatu * 180. / pi)
        ierr = nf90_put_var(nid, vid_latv, rlatv * 180. / pi)
        ierr = nf90_put_var(nid, vid_lev, presnivs)
        ierr = nf90_put_var(nid, vid_cu, cu)
        ierr = nf90_put_var(nid, vid_cv, cv)
        ierr = nf90_put_var(nid, vid_au, zu)
        ierr = nf90_put_var(nid, vid_av, zv)
        call nf95_put_var(nid, varid_alpha_t, zt)
        call nf95_put_var(nid, varid_alpha_q, zq)
! --------------------------------------------------------------------
! Cr�ation des variables sauvegard�es
! --------------------------------------------------------------------
        ierr = nf90_redef(nid)
! Pressure (GCM)
        dim4=(/id_lonv,id_latu,id_lev,id_tim/)
        ierr = nf90_def_var(nid,"SP",nf90_float,dim4,varid)
! Surface pressure (guidage)
        IF (guide_P) THEN
            dim3=(/id_lonv,id_latu,id_tim/)
            ierr = nf90_def_var(nid,"ps",nf90_float,dim3,varid)
        ENDIF
! Zonal wind
        IF (guide_u) THEN
            dim4=(/id_lonu,id_latu,id_lev,id_tim/)
            ierr = nf90_def_var(nid,"u",nf90_float,dim4,varid)
            ierr = nf90_def_var(nid,"ua",nf90_float,dim4,varid)
            ierr = nf90_def_var(nid,"ucov",nf90_float,dim4,varid)
        ENDIF
! Merid. wind
        IF (guide_v) THEN
            dim4=(/id_lonv,id_latv,id_lev,id_tim/)
            ierr = nf90_def_var(nid,"v",nf90_float,dim4,varid)
            ierr = nf90_def_var(nid,"va",nf90_float,dim4,varid)
            ierr = nf90_def_var(nid,"vcov",nf90_float,dim4,varid)
        ENDIF
! Pot. Temperature
        IF (guide_T) THEN
            dim4=(/id_lonv,id_latu,id_lev,id_tim/)
            ierr = nf90_def_var(nid,"teta",nf90_float,dim4,varid)
        ENDIF
! Specific Humidity
        IF (guide_Q) THEN
            dim4=(/id_lonv,id_latu,id_lev,id_tim/)
            ierr = nf90_def_var(nid,"q",nf90_float,dim4,varid)
        ENDIF

        ierr = nf90_enddef(nid)
        ierr = nf90_close(nid)
    ENDIF ! timestep=0

! --------------------------------------------------------------------
! Enregistrement du champ
! --------------------------------------------------------------------

    ierr=nf90_open("guide_ins.nc",nf90_write,nid)

    IF (varname=="SP") timestep=timestep+1

    ierr = nf90_inq_varid(nid,varname,varid)
    SELECT CASE (varname)
    CASE ("SP","ps")
        start=(/1,1,1,timestep/)
        count=(/iip1,jjp1,llm,1/)
    CASE ("v","va","vcov")
        start=(/1,1,1,timestep/)
        count=(/iip1,jjm,llm,1/)
    CASE DEFAULT
        start=(/1,1,1,timestep/)
        count=(/iip1,jjp1,llm,1/)
    END SELECT

!$OMP END MASTER
!$OMP BARRIER

    SELECT CASE (varname)

    CASE("u","ua")
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
            field_glo(:,2:jjm,l)=field_glo(:,2:jjm,l)/cu(:,2:jjm)
            field_glo(:,1,l)=0. ; field_glo(:,jjp1,l)=0.
        ENDDO
    CASE("v","va")
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
           field_glo(:,:,l)=field_glo(:,:,l)/cv(:,:)
        ENDDO
    END SELECT

!    if (varname=="ua") then
!    call dump2d(iip1,jjp1,field_glo,'ua gui1 1ere couche ')
!    call dump2d(iip1,jjp1,field_glo(:,:,llm),'ua gui1 llm ')
!    endif

!$OMP MASTER

    ierr = nf90_put_var(nid, varid, field_glo, start, count)
    ierr = nf90_close(nid)

       DEALLOCATE(field_glo)
!$OMP END MASTER
!$OMP BARRIER

  END SUBROUTINE guide_out


!===========================================================================
  subroutine correctbid(iim,nl,x)
    integer iim,nl
    real x(iim+1,nl)
    integer i,l
    real zz

    do l=1,nl
        do i=2,iim-1
            if(abs(x(i,l)).gt.1.e10) then
               zz=0.5*(x(i-1,l)+x(i+1,l))
              print*,'correction ',i,l,x(i,l),zz
               x(i,l)=zz
            endif
         enddo
     enddo
     return
  end subroutine correctbid


!====================================================================
! Ascii debug output. Could be reactivated
!====================================================================

subroutine dump2du(var,varname)
     use parallel_lmdz
use mod_hallo
USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
implicit none



      CHARACTER (len=*) :: varname


real, dimension(ijb_u:ije_u) :: var

real, dimension(ip1jmp1) :: var_glob

    RETURN

    call barrier
    CALL Gather_field_u(var,var_glob,1)
    call barrier

    if (mpi_rank==0) then
       call dump2d(iip1,jjp1,var_glob,varname)
    endif

    call barrier

    return
    end subroutine dump2du

!====================================================================
! Ascii debug output. Could be reactivated
!====================================================================
subroutine dumpall
     USE dimensions_mod, ONLY: iim, jjm, llm, ndm
     USE paramet_mod_h
     USE comgeom_mod_h
implicit none


     call barrier
     call dump2du(alpha_u(ijb_u:ije_u),'  alpha_u couche 1')
     call dump2du(unat2(:,jjbu:jjeu,nlevnc),'  unat2 couche nlevnc')
     call dump2du(ugui1(ijb_u:ije_u,1)*sqrt(unscu2(ijb_u:ije_u)),'  ugui1 couche 1')
     return
end subroutine dumpall

!===========================================================================
END MODULE guide_loc_mod
