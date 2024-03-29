!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/CCTM/src/vdiff/acm2/vdiffacm2.F,v 1.13 2012/01/19 14:37:47 yoj Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VDIFF ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Asymmetric Convective Model v2 (ACM2) -- Pleim(2006)
C Function:
C   calculates and writes dry deposition.
C   calculates vertical diffusion

C Subroutines and Functions Called:
C   INIT3, SEC2TIME, TIME2SEC, WRITE3, NEXTIME,
C   M3EXIT, EDDYX, TRI, MATRIX, PA_UPDATE_EMIS, PA_UPDATE_DDEP

C Revision History:
C   Analogous to VDIFFIM (Eddy diffusion PBL scheme)

C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module (GRID_CONF)
C    7 Jun 05 P.Bhave: added call to OPSSEMIS if MECHNAME='AE4';
C                      added TSTEP to RDEMIS_AE call vector
C    Aug 05 J. Pleim Update to v4.5
C
C    Jan 06 J. Pleim ACM2 implementation
C    Jul 06 J.Young: emis module
C    Nov 06 B.Hutzell, J.Young: depv module
C    Mar 08 G.Sarwar: added heterogeneous reaction for HONO. It affects HONO, NO2 and HNO3
C    Apr 08 J.Young: optional inline beis and plume rise emissions
C 21 Jun 10 J.Young: convert for Namelist redesign
C 16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C 10 May 11 J.Bash: updated for NH3 bidirectional exchange
C 11 May 11 D.Wong: - incorporated twoway model implementation
C 20 May 11 D.Schwede: modified to include mosaic
C 15 Aug 11 D.Schwede and G. Sarwar: Modified HONO calculation so that deposition
C           velocity for NO2 that is output in DEPV file does not include the loss due
C           to the heterogeneous reaction.
C           This additional loss is now accounted for in vdiff.F
C 14 Sep 11 J.Young: fix TSTEP bug in call to DEPV_INIT
C 11 Apr 13 J.Young: fix double adjustment of conc for DDBF in heterogeneous HONO
C           if-then-else clauses; eliminate some white space
C 13 May 13 J.Young: access met data from VDIFF_MET module
C                    change CRANKP to THBAR, CRANKQ to THETA
C 25 May 13 J.Young: re-do the acm/eddy algorithm for computational efficiency
C 30 Apr 14 J.Young: switch THBAR and THETA
C  2 May 14 J.Pleim, J.Young: replace the banded tridiagonal matrix solver for the
C           convective PBL, with the ACM1 matrix solver followed by the tridiagonal
C           matrix solver
C 30 May 14 J.Young: split vdiff calculation out of vdiff proc.
C 31 Jul 14 J.Young: incorporate gravitational settling for J-,K-mode aerosols
C-----------------------------------------------------------------------

      USE CGRID_SPCS          ! CGRID mechanism species
      USE GRID_CONF
      USE EMIS_DEFN          
      USE DEPV_DEFN
      USE VDIFF_MET
      USE VDIFF_MAP
      USE UTILIO_DEFN
      USE BIDI_MOD
      USE LSM_MOD, ONLY: n_lufrac
      USE VDIFF_DIAG
#ifdef parallel
      USE SE_MODULES            ! stenex (using SE_GLOBAL_SUM_MODULE)
#else
      USE NOOP_MODULES          ! stenex (using NOOP_GLOBAL_SUM_MODULE)
#endif

      IMPLICIT NONE

      INCLUDE SUBST_PACTL_ID  ! PA control parameters
      INCLUDE SUBST_CONST     ! constants
      INCLUDE SUBST_FILES_ID  ! file name parameters

      CHARACTER( 120 ) :: XMSG = ' '

C Arguments:

      REAL, POINTER :: CGRID( :,:,:,: )              !  concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 3 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)
                                ! TSTEP(3) = twoway model time step w.r.t. wrf time
                                !            step and wrf/cmaq call frequency

C Parameters:

C External Functions: None

C Local Variables: 
      CHARACTER( 16 ), SAVE :: PNAME = 'VDIFFPROC'
      CHARACTER( 16 ), SAVE :: AERO_GRAV_SETL = 'CTM_GRAV_SETL'
      CHARACTER( 80 ) :: VARDESC                ! env variable description
      LOGICAL, SAVE :: GRAV_SETL
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: WSTEP  = 0               ! local write counter
      INTEGER  STATUS                           ! ENV... status

      REAL          :: FCJACMF( NCOLS,NROWS,NLAYS )  ! 1/ mid-full layer vert Jac factor
      REAL             LRDX3M                        ! loop local RDX3M( L )
      REAL             FCMSF                         ! loop local RMSFX4( C,R )

      REAL, ALLOCATABLE, SAVE :: CNGRD( :,:,:,: )     ! cgrid replacement
      REAL, ALLOCATABLE, SAVE :: DDEP     ( :,:,: )   ! ddep accumulator
      REAL, ALLOCATABLE, SAVE :: ICMP     ( :,:,: )   ! component flux accumlator 
      REAL, ALLOCATABLE, SAVE :: DDEPJ    ( :,:,:,: ) ! ddep for mosaic
      REAL, ALLOCATABLE, SAVE :: DDEPJ_FST( :,:,:,: ) ! ddep for stomtal/cuticular pathway

      REAL     :: WRDD( NCOLS,NROWS )                 ! ddep write buffer
      REAL     :: WRDDJ( NCOLS,NROWS,N_LUFRAC+1 )     ! mosaic ddep write buffer
      REAL     :: WRDDJ_FST( NCOLS,NROWS,N_LUFRAC+1 ) ! mosaic stomatal flux write buffer

      REAL, ALLOCATABLE, SAVE :: DDEP_PA  ( :,:,: )   ! ddep for process analysis
      REAL, ALLOCATABLE, SAVE :: EMIS_PA( :,:,:,: )   ! emis for process analysis

      INTEGER, SAVE :: N_SPC_CGRID              ! no. of CGRID species

      REAL     :: EDDYV ( NCOLS,NROWS,NLAYS )   ! from EDYINTB
      REAL     :: SEDDY ( NLAYS,NCOLS,NROWS )   ! flipped EDDYV
      REAL        DTSEC                         ! model time step in seconds

      REAL, ALLOCATABLE, SAVE :: VSED_AE( :,:,:,: )

C Local Variables

      INTEGER, SAVE :: LOGDEV
 
      INTEGER     ASTAT
      INTEGER     C, R, L, S, V, I, J, OFF      ! loop induction variables
      INTEGER     MDATE, MTIME, MSTEP           ! internal simulation date&time

      INTERFACE
         SUBROUTINE PA_UPDATE_EMIS ( PNAME, VDEMIS, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN )  :: PNAME
            REAL,           INTENT( IN )  :: VDEMIS( :,:,:,: )
            INTEGER,        INTENT( IN )  :: JDATE, JTIME
            INTEGER,        INTENT( IN )  :: TSTEP( 3 )
         END SUBROUTINE PA_UPDATE_EMIS
         SUBROUTINE PA_UPDATE_DDEP ( PNAME, DDEP, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN )  :: PNAME
            REAL,           INTENT( IN )  :: DDEP( :,:,: )
            INTEGER,        INTENT( IN )  :: JDATE, JTIME
            INTEGER,        INTENT( IN )  :: TSTEP( 3 )
         END SUBROUTINE PA_UPDATE_DDEP
         SUBROUTINE CONV_CGRID ( CGRID, JDATE, JTIME, CNGRD )
            IMPLICIT NONE
            REAL, POINTER :: CGRID( :,:,:,: )
            INTEGER,        INTENT( IN )  :: JDATE, JTIME
            REAL,           INTENT( INOUT ) :: CNGRD( :,:,:,: )
         END SUBROUTINE CONV_CGRID
         SUBROUTINE REV_CGRID ( CNGRD, JDATE, JTIME, CGRID )
            IMPLICIT NONE
            REAL,           INTENT( INOUT ) :: CNGRD( :,:,:,: )
            INTEGER,        INTENT( IN )  :: JDATE, JTIME
            REAL, POINTER :: CGRID( :,:,:,: )
         END SUBROUTINE REV_CGRID
         SUBROUTINE EDDYX ( EDDYV )
            USE GRID_CONF
            USE VDIFF_MET
            USE UTILIO_DEFN
#ifdef parallel
            USE SE_MODULES
#else
            USE NOOP_MODULES
#endif
            IMPLICIT NONE
            REAL,           INTENT( OUT ) :: EDDYV( :,:,: )
         END SUBROUTINE EDDYX
         SUBROUTINE VDIFFACMX( DTSEC, SEDDY, DDEP, ICMP, DDEPJ, DDEPJ_FST, CNGRD )
            USE CGRID_SPCS
            USE GRID_CONF
            USE EMIS_DEFN
            USE DEPV_DEFN
            USE VDIFF_MET
            USE VDIFF_MAP
            USE UTILIO_DEFN
            USE BIDI_MOD
            USE LSM_MOD, ONLY: n_lufrac
            IMPLICIT NONE
            REAL, INTENT( IN )    :: DTSEC
            REAL, INTENT( INOUT ) :: SEDDY    ( :,:,: )
            REAL, INTENT( INOUT ) :: DDEP     ( :,:,: )
            REAL, INTENT( INOUT ) :: ICMP     ( :,:,: )
            REAL, INTENT( INOUT ) :: DDEPJ    ( :,:,:,: )
            REAL, INTENT( INOUT ) :: DDEPJ_FST( :,:,:,: )
            REAL, INTENT( INOUT ) :: CNGRD    ( :,:,:,: )
         END SUBROUTINE VDIFFACMX
         SUBROUTINE SEDI( DTSEC, SEDVEL, CGRID, CNGRD )
            USE CGRID_SPCS
            USE GRID_CONF
            USE VDIFF_MET
            USE VDIFF_MAP
            USE UTILIO_DEFN
            IMPLICIT NONE
            REAL, INTENT( IN )    :: DTSEC
            REAL, INTENT( OUT )   :: SEDVEL   ( :,:,:,: )
            REAL, POINTER         :: CGRID    ( :,:,:,: )
            REAL, INTENT( INOUT ) :: CNGRD    ( :,:,:,: )
         END SUBROUTINE SEDI

      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

         IF ( .NOT. DEPV_INIT ( JDATE, JTIME, TSTEP, CGRID ) ) THEN
            XMSG = 'Failure initializing deposition velocities module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C create global maps
         IF ( .NOT. VDIFF_MAP_INIT( N_SPC_DEPV ) ) THEN
            XMSG = 'Failure initializing index mapping module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Initialize the met data
         IF ( .NOT. VDIFF_MET_INIT( JDATE, JTIME ) ) THEN
            XMSG = 'Failure initializing met module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get gravitational settling (sedi) flag.
         GRAV_SETL = .TRUE.         ! default
         VARDESC = 'Using J-,K-mode aerosols gravitational settling'
         GRAV_SETL = ENVYN( AERO_GRAV_SETL, VARDESC, GRAV_SETL, STATUS )
         IF ( STATUS .EQ. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC

C Get diagnostic files flag.
         VDIFFDIAG = .FALSE.         ! default
         VARDESC = 'Writing the VDIFF diagnostic files'
         VDIFFDIAG = ENVYN( VDIFF_DIAG_FILE, VARDESC, VDIFFDIAG, STATUS )
         IF ( STATUS .EQ. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
 
C Set output file characteristics based on COORD.EXT and open the dry dep file
         IF ( MYPE .EQ. 0 ) CALL OPDDEP ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP, ABFLUX )
         
         IF ( ABFLUX .AND. MYPE .EQ. 0 ) CALL OPASX_MEDIA( JDATE, JTIME, TSTEP( 1 ) )

C Allocate and initialize dry deposition array

         ALLOCATE ( DDEP( N_SPC_DEPV,NCOLS,NROWS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DDEP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         DDEP = 0.0   ! array assignment
  
         ALLOCATE ( ICMP( LCMP,NCOLS,NROWS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ICMP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         ICMP = 0.0   ! array assignment

         IF ( .NOT. EMIS_INIT ( JDATE, JTIME, TSTEP( 1 ) ) ) THEN
            XMSG = 'Failure initializing emissions module'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Set up for process analysis
         IF ( LIPR ) THEN
            ALLOCATE ( EMIS_PA( NCOLS,NROWS,EMLAYS,N_SPC_EMIS+1 ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'EMIS_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            ALLOCATE ( DDEP_PA( NCOLS,NROWS,N_SPC_DEPV ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'DDEP_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

C Set up for grav. settling
         IF ( GRAV_SETL ) THEN
            ALLOCATE ( VSED_AE( N_AE_SPC,NLAYS,NCOLS,NROWS ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating VSED_AE'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

C Open vdiff diagnostics file (ioapi header from cgrd)
         IF ( VDIFFDIAG ) THEN
            IF ( .NOT. VDIFF_DIAG_INIT ( JDATE, JTIME, TSTEP( 1 ), GRAV_SETL ) ) THEN
               XMSG = 'Failure initializing vdiff diagnostics module'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         N_SPC_CGRID = SIZE ( CGRID,4 )

         ALLOCATE ( CNGRD( N_SPC_CGRID,NLAYS,NCOLS,NROWS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CNGRD'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         CNGRD = 0.0   ! array assignment
 
         IF ( MOSAIC ) THEN
            ALLOCATE ( DDEPJ( N_LUFRAC,N_SPC_DEPV,NCOLS,NROWS ), STAT = ASTAT )
            IF ( ASTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DDEPJ'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            DDEPJ = 0.0   ! array assignment
            IF ( MYPE .EQ. 0 ) CALL OPDDEP_MOS ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP )
            IF ( FST ) THEN
               ALLOCATE ( DDEPJ_FST( N_LUFRAC,N_SPC_DEPV,NCOLS,NROWS ), STAT = ASTAT )
               IF ( ASTAT .NE. 0 ) THEN
                  XMSG = 'Failure allocating DDEPJ_FST'
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
               DDEPJ_FST = 0.0   ! array assignment
               IF ( MYPE .EQ. 0 ) CALL OPDDEP_FST ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP )
            END IF   ! if Fst
         END IF   ! if Mosaic

      END IF   !  if Firstime

      MDATE = JDATE
      MTIME = JTIME
      MSTEP = TIME2SEC( TSTEP( 2 ) )
      DTSEC = FLOAT( MSTEP )
      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

C Convert non-molar mixing ratio species and re-order CGRID
      CALL CONV_CGRID ( CGRID, MDATE, MTIME, CNGRD )

C read & interpolate deposition velocities
      CALL GET_DEPV ( MDATE, MTIME, TSTEP, CGRID )
      
C read & interpolate met data
      CALL GET_MET ( MDATE, MTIME, MSTEP )

      IF ( GRAV_SETL ) THEN
C Get gravitational settling velocity for the vsed aero species:
C AERO_SEDV assumes that every aero species is dry deposited and is diffused (trns)
C Calculate the changes in the layer J-,K-mode aerosol concentrations
         CALL SEDI( DTSEC, VSED_AE, CGRID, CNGRD )
      END IF

C read & interpolate emissions data => VDEMIS from EMIS_DEFN module
      CALL GET_EMIS ( MDATE, MTIME, TSTEP, CONVPA, RDX3F, RRHOJ, RJACM )

      IF ( LIPR ) THEN
         DO S = 1, N_SPC_EMIS+1
            DO L = 1, EMLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     EMIS_PA( C,R,L,S ) = VDEMIS( S,L,C,R )
                  END DO
               END DO
            END DO
         END DO
         CALL PA_UPDATE_EMIS ( 'VDIF', EMIS_PA, JDATE, JTIME, TSTEP )
      END IF

      CALL EDDYX ( EDDYV )

C EDDYV returned = Kz, where Kz is in m**2/sec

      DO L = 1, NLAYS
         LRDX3M = RDX3M( L )
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               FCJACMF( C,R,L ) = LRDX3M * RJACM( C,R,L ) * RJACF( C,R,L )
            END DO
         END DO
      END DO
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            FCMSF = RMSFX4( C,R )
            DO L = 1, NLAYS
               SEDDY( L,C,R ) = FCMSF * FCJACMF( C,R,L ) * EDDYV( C,R,L )
            END DO
         END DO
      END DO

      IF ( WSTEP .EQ. 0 ) THEN
         DDEP = 0.0                      ! array assignment
         ICMP = 0.0                      ! array assignment
         IF ( MOSAIC ) THEN
            DDEPJ = 0.0                  ! array assignment
            IF ( FST ) DDEPJ_FST = 0.0   ! array assignment
         END IF
      END IF

C Calculate the change in concentration and dry dep from vertical diffusion and vsed
      CALL VDIFFACMX( DTSEC, SEDDY, DDEP, ICMP, DDEPJ, DDEPJ_FST, CNGRD )

      IF ( VDIFFDIAG ) THEN
         NTICS = NTICS + 1
         NLPCR_SUM = NLPCR_SUM + NLPCR_MEAN    ! array assignment
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               NLPCR_MAX( C,R ) = MAX( NLPCR_MEAN( C,R ), NLPCR_MAX( C,R ) )
               NLPCR_MIN( C,R ) = MIN( NLPCR_MEAN( C,R ), NLPCR_MIN( C,R ) )
            END DO
         END DO
         IF ( GRAV_SETL ) THEN
            DTCCR_SUM = DTCCR_SUM + DTCCR_MEAN    ! array assignment
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  DTCCR_MAX( C,R ) = MAX( DTCCR_MEAN( C,R ), DTCCR_MAX( C,R ) )
                  DTCCR_MIN( C,R ) = MIN( DTCCR_MEAN( C,R ), DTCCR_MIN( C,R ) )
               END DO
            END DO
         END IF
      END IF

C Revert non-molar mixing ratio species and re-order CGRID
      CALL REV_CGRID ( CNGRD, MDATE, MTIME, CGRID )

C If last call this hour:  write accumulated depositions:

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      IF ( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) THEN
         MDATE = JDATE
         MTIME = JTIME
         CALL NEXTIME( MDATE, MTIME, TSTEP( 2 ) )
         WSTEP = 0

         DO V = 1, N_SPC_DDEP
            S = DD2DV( V )
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  WRDD( C,R ) = DDEP( S,C,R )
               END DO
            END DO

            IF ( .NOT. WRITE3( CTM_DRY_DEP_1, DDEP_SPC( V ),
     &                 MDATE, MTIME, WRDD ) ) THEN
               XMSG = 'Could not write ' // CTM_DRY_DEP_1 // ' file'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
    
            IF ( ABFLUX .AND. TRIM( DDEP_SPC( V ) ) .EQ. 'NH3' ) THEN
               DO I = 1, LCMP 
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                        WRDD( C,R ) = ICMP( I,C,R )
                     END DO
                  END DO
                  IF ( .NOT. WRITE3( CTM_DRY_DEP_1, CMPSPC( I ),
     &                 MDATE, MTIME, WRDD ) ) THEN
                     XMSG = 'Could not write ' // CTM_DRY_DEP_1 // ' file'
                     CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                  END IF     
               END DO       
            ENDIF

         END DO 

         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CTM_DRY_DEP_1,
     &         'for date and time', MDATE, MTIME
         
C Write vdiff diagnostics
         IF ( VDIFFDIAG ) THEN
            IF ( GRAV_SETL ) THEN   ! Write vsed diagnostics

               DO V = 1, N_VSED
                  S = VSED_MAP( V )
                  DO L = 1, NLAYS
                     DO R = 1, MY_NROWS
                        DO C = 1, MY_NCOLS
                           VSED_BUF( C,R,L,V ) = VSED_AE( S,L,C,R )
                        END DO
                     END DO
                  END DO
                  IF ( .NOT. WRITE3( CTM_VSED_DIAG, VSED_NAME( V ),
     &                               MDATE, MTIME, VSED_BUF( 1,1,1,V ) ) ) THEN
                     XMSG = 'Could not write ' // TRIM( VSED_NAME( V ) )
     &                    // ' to ' // CTM_VSED_DIAG
                     CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                  END IF
               END DO

               WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &               'Timestep written to', CTM_VSED_DIAG,
     &               'for date and time', MDATE, MTIME

            END IF   ! GRAV_SETL

C Write other diagnostics
            NLPCR_MEAN = NLPCR_SUM / FLOAT( NTICS )
            IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'NLP_MEAN',
     &                         MDATE, MTIME, NLPCR_MEAN ) ) THEN
               XMSG = 'Could not write ' //  'NLP_MEAN to ' // CTM_VDIFF_DIAG
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'NLP_MAX',
     &                         MDATE, MTIME, NLPCR_MAX ) ) THEN
               XMSG = 'Could not write ' //  'NLP_MAX to ' // CTM_VDIFF_DIAG
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'NLP_MIN',
     &                         MDATE, MTIME, NLPCR_MIN ) ) THEN
               XMSG = 'Could not write ' //  'NLP_MIN to ' // CTM_VDIFF_DIAG
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            NLPCR_MAX = 0.0      ! array assignment
            NLPCR_MIN = 9.9E30   ! array assignment
            NLPCR_SUM = 0.0      ! array assignment

            IF ( GRAV_SETL ) THEN   ! Write vsed diagnostics
               DTCCR_MEAN = DTCCR_SUM / FLOAT( NTICS )
               IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'SEDI_DTC_MEAN',
     &                            MDATE, MTIME, DTCCR_MEAN ) ) THEN
                  XMSG = 'Could not write ' //  'SEDI_DTC_MEAN to ' // CTM_VDIFF_DIAG
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
               IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'SEDI_DTC_MAX',
     &                            MDATE, MTIME, DTCCR_MAX ) ) THEN
                  XMSG = 'Could not write ' //  'SEDI_DTC_MAX to ' // CTM_VDIFF_DIAG
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
               IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'SEDI_DTC_MIN',
     &                            MDATE, MTIME, DTCCR_MIN ) ) THEN
                  XMSG = 'Could not write ' //  'SEDI_DTC_MIN to ' // CTM_VDIFF_DIAG
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF
               DTCCR_MAX = 0.0      ! array assignment
               DTCCR_MIN = 9.9E30   ! array assignment
               DTCCR_SUM = 0.0      ! array assignment
            END IF

            CNVCT = 0.0   ! array assignment
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  IF ( CONVCT( C,R ) ) CNVCT( C,R ) = 1.0
               END DO
            END DO
            IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'CONVCT',
     &                         MDATE, MTIME, CNVCT ) ) THEN
               XMSG = 'Could not write ' //  'convct to ' // CTM_VDIFF_DIAG
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

            IF ( .NOT. WRITE3( CTM_VDIFF_DIAG, 'LPBL',
     &                         MDATE, MTIME, REAL( LPBL ) ) ) THEN
               XMSG = 'Could not write ' //  'lpbl to ' // CTM_VDIFF_DIAG
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

            WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6, I6 )' )
     &            'Timestep written to', CTM_VDIFF_DIAG,
     &            'for date and time (and ntics)', MDATE, MTIME, NTICS
            NTICS = 0

         END IF
         
         IF ( MOSAIC ) THEN

            DO V = 1, N_SPC_DDEP
               S = DD2DV( V )
               WRDD = 0.0 ! reuse array since it has already been written for hour
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     DO J = 1, N_LUFRAC
                        WRDD( C,R ) = WRDD( C,R ) + DDEPJ( J,S,C,R ) * LUFRAC( J,C,R )
                        WRDDJ( C,R,J ) = DDEPJ( J,S,C,R )
                     END DO
                     WRDDJ( C,R,N_LUFRAC+1 ) = WRDD( C,R )  ! last array element is total across all land use categories
                  END DO
               END DO

               IF ( .NOT. WRITE3( CTM_DRY_DEP_MOS, DDEP_SPC( V ),
     &                     MDATE, MTIME, WRDDJ ) ) THEN
                  XMSG = 'Could not write ' // CTM_DRY_DEP_MOS // ' file'
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF

            END DO

            WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &             'Timestep written to', CTM_DRY_DEP_MOS,
     &             'for date and time', MDATE, MTIME

            IF ( FST ) THEN

               DO V = 1, N_SPC_DDEP
                  S = DD2DV( V )
                  WRDD = 0.0 ! reuse array since it has already been written for hour
                  DO R = 1, MY_NROWS
                     DO C = 1, MY_NCOLS
                        DO J = 1, N_LUFRAC
                           WRDD( C,R ) = WRDD( C,R ) + DDEPJ_FST( J,S,C,R ) * LUFRAC( J,C,R )
                           WRDDJ_FST( C,R,J ) = DDEPJ_FST( J,S,C,R )
                           IF ( DDEPJ_FST( J,S,C,R ) .GT. DDEPJ( J,S,C,R ) ) THEN
                              WRITE( LOGDEV,* ) 'FST too big !!!'
                              WRITE( LOGDEV,* ) 'J,S,C,R = ', J, S, C, R
                              WRITE( LOGDEV,* ) 'DDEPJ,DDEPJ_FST: ', DDEPJ( J,S,C,R ), DDEPJ_FST( J,S,C,R )
                              WRITE( LOGDEV,* ) 'DDEP Species: ', DDEP_SPC( V )
                              WRITE( LOGDEV,* ) 'Time and date: ', MTIME, MDATE
                           END IF
                        END DO
                        WRDDJ_FST( C,R,N_LUFRAC+1 ) = WRDD( C,R )  ! last array element is total across all land use categories
                    END DO
                  END DO
             
                  IF ( .NOT. WRITE3( CTM_DRY_DEP_FST, DDEP_SPC( V ),
     &                       MDATE, MTIME, WRDDJ_FST ) ) THEN
                     XMSG = 'Could not write ' // CTM_DRY_DEP_FST // ' file'
                     CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                  END IF
             
               END DO

               WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &               'Timestep written to', CTM_DRY_DEP_FST,
     &               'for date and time', MDATE, MTIME

            END IF   ! FST

         END IF   ! MOSAIC

         IF ( ABFLUX ) THEN    
            CALL WRASX_MEDIA( MDATE, MTIME )
         END IF 
 
         IF ( LIPR ) THEN
            DO V = 1, N_SPC_DEPV
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     DDEP_PA( C,R,V ) = DDEP( V,C,R )
                  END DO
               END DO
            END DO
            CALL PA_UPDATE_DDEP ( 'VDIF', DDEP_PA, JDATE, JTIME, TSTEP )
         END IF

C re-set dry deposition array to zero

         DDEP = 0.0
         ICMP = 0.0
         IF ( MOSAIC ) THEN
            DDEPJ = 0.0   ! array assignment
            IF ( FST ) DDEPJ_FST = 0.0   ! array assignment
         END IF

      END IF

      RETURN
      END
