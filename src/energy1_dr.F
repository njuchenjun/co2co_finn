! --------------------------------------------------------------------
!     input  Q:    CARTESIAN IN ANGSTROM
!     return VOUT: ENERGY IN eV
!     return DVDQ: DERIVATIVES IN eV/Angstrom
! --------------------------------------------------------------------
#     DEFINE NATOM 5
#     DEFINE NBOND (NATOM*(NATOM-1)/2)
! --------------------------------------------------------------------
      SUBROUTINE ENERGY1(Q,VOUT,DVDQ)
      IMPLICIT NONE
      REAL*8,INTENT(IN) :: Q(3,NATOM)
      REAL*8,INTENT(OUT) :: VOUT,DVDQ(3,NATOM)

      INTEGER I,J,K
      REAL*8 :: R(NBOND),DR(3,NATOM,NBOND),DVDR(NBOND)
      INTEGER,SAVE :: MB(NBOND),NB(NBOND)
      INTEGER,SAVE :: INIT=0

      REAL*8,PARAMETER :: TOANG=0.5291772083D0 !1bohr=0.5291772083a ! from a.u.
      REAL*8,PARAMETER :: cm=219474.5d0     !1hartree=219474.5cm-1
      REAL*8,PARAMETER :: ev=27.21138505d0

      real*8 :: v,vco2,vco,vint
      real*8 :: dvco,dvco2(3)

      IF(INIT.EQ.0) THEN
        INIT=1
        K=1
        DO I=1,NATOM-1
        DO J=I+1,NATOM
          MB(K)=I
          NB(K)=J
          K=K+1
        ENDDO
        ENDDO
      ENDIF

      CALL INTERN(Q,MB,NB,R,DR)

      call pes_co2co(Q,v,vco2,vco,vint,dvdr,dvco2,dvco) !v : cm-1


      vout=(vint/cm + vco2 + vco/cm)*ev
      
      !dvdr: cm-1/angstrom, dvco2: a.u., dco: cm-1/bohr
      dvdr(1)=dvdr(1) + dvco2(3)*cm/TOANG !! now dvdr in unit of cm-1/Angstrom
      dvdr(2)=dvdr(2) + dvco2(1)*cm/TOANG
      dvdr(5)=dvdr(5) + dvco2(2)*cm/TOANG
      dvdr(10)=dvdr(10) + dvco/TOANG

      dvdr = dvdr/cm*ev ! now in eV/Angstrom

      CALL DVDX(MB,NB,DR,DVDR,DVDQ) !! return DVDQ in unit of eV/Angstrom

      RETURN
      END
! --------------------------------------------------------------------
      SUBROUTINE INTERN(Q,MB,NB,R,DR)
! ---- CALCULATES BONDLENGTHS AND THEIR DERIVATIVES
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8,INTENT(IN) :: Q(3,NATOM)
      INTEGER,INTENT(IN) :: MB(NBOND),NB(NBOND)
      REAL*8,INTENT(OUT) :: R(NBOND),DR(3,NATOM,NBOND)
! ---- BOND LENGTHS
      DO I=1,NBOND  !amount of bond length change
        R(I)=DSQRT(
     .   (Q(1,MB(I))-Q(1,NB(I)))*(Q(1,MB(I))-Q(1,NB(I)))
     .  +(Q(2,MB(I))-Q(2,NB(I)))*(Q(2,MB(I))-Q(2,NB(I)))
     .  +(Q(3,MB(I))-Q(3,NB(I)))*(Q(3,MB(I))-Q(3,NB(I))))
      ENDDO
! ---- BOND LENGTH DERIVATIVES
      DO I=1,NBOND  !x,y,z direction partial derivative about R
         DR(1,MB(I),I)=(Q(1,MB(I))-Q(1,NB(I)))/R(I)
         DR(2,MB(I),I)=(Q(2,MB(I))-Q(2,NB(I)))/R(I)
         DR(3,MB(I),I)=(Q(3,MB(I))-Q(3,NB(I)))/R(I)
         DR(1,NB(I),I)=-DR(1,MB(I),I)
         DR(2,NB(I),I)=-DR(2,MB(I),I)
         DR(3,NB(I),I)=-DR(3,MB(I),I)
      ENDDO
      RETURN
      END SUBROUTINE INTERN
! --------------------------------------------------------------------
      SUBROUTINE DVDX(MB,NB,DR,DVDR,DV)
! ---- CALCULATES DV/DX FROM DV/DR*DR/DX
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MB(NBOND),NB(NBOND)
      REAL*8,INTENT(IN) :: DR(3,NATOM,NBOND),DVDR(NBOND)
      REAL*8,INTENT(OUT) :: DV(3,NATOM)
      INTEGER :: I
      DV=0.D0
      DO I=1,NBOND
         DV(1,MB(I))=DV(1,MB(I))+DVDR(I)*DR(1,MB(I),I)
         DV(1,NB(I))=DV(1,NB(I))+DVDR(I)*DR(1,NB(I),I)
         DV(2,MB(I))=DV(2,MB(I))+DVDR(I)*DR(2,MB(I),I)
         DV(2,NB(I))=DV(2,NB(I))+DVDR(I)*DR(2,NB(I),I)
         DV(3,MB(I))=DV(3,MB(I))+DVDR(I)*DR(3,MB(I),I)
         DV(3,NB(I))=DV(3,NB(I))+DVDR(I)*DR(3,NB(I),I)
      ENDDO
      RETURN
      END SUBROUTINE DVDX
! --------------------------------------------------------------------
