!=======================================================================
! ABAQUS/STANDARD USER DEFINED FIELD SUBROUTINE: USDFLD
!=======================================================================
! Project:    05. Secondary Transmission Shaft Straightening & Springback
! Model CAD:  husson_shaft.step (Z = -8.6 to 285.5 mm, L = 294.1 mm)
! Subroutine: USDFLD (User Defined Field)
! Solver:     Abaqus/Standard (Implicit)
! Language:   Fortran 2008 / Fixed-Free Form Compatible
! Author:     Victor Maia (Simulation & FEA Specialist)
!
! Physical Formulation:
! Evaluates through-thickness radial depth from the outer surface of the 
! stepped secondary transmission shaft:
!     d(r, z) = R_outer(z) - sqrt((x - x0)^2 + (y - y0)^2)
!
! Outer radius profile R_outer(z) is mapped according to the CAD 
! geometry of the secondary shaft (husson_shaft.step):
!   - Bearing Seat 1 (Z < 35 mm):          R = 12.5 mm  (Dia 25 mm)
!   - Main Central Span (35 <= Z <= 180):   R = 19.0 mm  (Dia 38 mm) [Punch zone]
!   - Stepped Shoulder (180 < Z <= 250):   R = 25.5 mm  (Dia 51 mm)
!   - Bearing Seat 2 (Z > 250 mm):         R = 15.0 mm  (Dia 30 mm)
!
! Maps depth d(r, z) directly to Field Variable 1 (FIELD(1)) at every 
! material integration point. Abaqus/Standard continuously interpolates 
! the depth-dependent Chaboche nonlinear kinematic hardening parameters
! and yield strength through *PLASTIC, DEPENDENCIES=1.
!
! State Variables (DEPVAR = 3):
!   STATEV(1) = Radial depth from outer case surface [mm]
!   STATEV(2) = Equivalent plastic strain (PEEQ) [-]
!   STATEV(3) = Current radial distance from centerline [mm]
!=======================================================================

      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     & TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     & KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)

      INCLUDE 'ABA_PARAM.INC'

      ! Argument Declarations
      CHARACTER*80 CMNAME, ORNAME
      CHARACTER*3  FLGRAY(15)
      DOUBLE PRECISION FIELD(NFIELD)
      DOUBLE PRECISION STATEV(NSTATV)
      DOUBLE PRECISION DIRECT(3,3), T(3,3), TIME(2), COORD(3)
      DOUBLE PRECISION PNEWDT, CELENT, DTIME
      INTEGER NFIELD, NSTATV, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      INTEGER NDI, NSHR, LACCFLA
      INTEGER JMAC(*), JMATYP(*)

      ! Local Variables & Arrays for GETVRM
      DOUBLE PRECISION ARRAY(15)
      INTEGER JARRAY(15)
      INTEGER JRCD

      ! Shaft Centerline Coordinates (Global Z-axis alignment)
      DOUBLE PRECISION, PARAMETER :: X0 = 0.0D0
      DOUBLE PRECISION, PARAMETER :: Y0 = 0.0D0

      DOUBLE PRECISION :: X_CURR, Y_CURR, Z_CURR
      DOUBLE PRECISION :: R_CURR, R_OUTER, DEPTH, PEEQ

      !-----------------------------------------------------------------
      ! 1. Query Current Material Integration Point Coordinates
      !-----------------------------------------------------------------
      CALL GETVRM('COORD', ARRAY, JARRAY, FLGRAY, JRCD, JMAC, JMATYP,
     &            MATLAYO, LACCFLA)
      
      IF (JRCD .NE. 0) THEN
         WRITE(6, 100) NOEL, NPT, KSTEP, KINC
 100     FORMAT(' *** ERROR in USDFLD: GETVRM failed to retrieve COORD',
     &          ' for Element ', I8, ' at Int. Point ', I3, 
     &          ' (Step ', I3, ', Inc ', I6, ')')
         CALL XIT
      END IF

      X_CURR = ARRAY(1)
      Y_CURR = ARRAY(2)
      Z_CURR = ARRAY(3)

      !-----------------------------------------------------------------
      ! 2. Piecewise Outer Radius Mapping (Stepped Shaft Profile)
      !    Calibrated to husson_shaft.step
      !-----------------------------------------------------------------
      IF (Z_CURR .LT. 35.0D0) THEN
         R_OUTER = 12.5D0
      ELSE IF (Z_CURR .LE. 180.0D0) THEN
         ! Main central span where press straightening is applied
         R_OUTER = 19.0D0
      ELSE IF (Z_CURR .LE. 250.0D0) THEN
         R_OUTER = 25.5D0
      ELSE
         R_OUTER = 15.0D0
      END IF

      !-----------------------------------------------------------------
      ! 3. Compute Radial Distance and Depth from Case Surface
      !-----------------------------------------------------------------
      R_CURR = SQRT((X_CURR - X0)**2 + (Y_CURR - Y0)**2)
      DEPTH  = R_OUTER - R_CURR

      ! Numerical boundary safeguard against mesh boundary tolerances
      IF (DEPTH .LT. 0.0D0) THEN
         DEPTH = 0.0D0
      END IF

      !-----------------------------------------------------------------
      ! 4. Query Equivalent Plastic Strain (PEEQ) History
      !-----------------------------------------------------------------
      CALL GETVRM('PEEQ', ARRAY, JARRAY, FLGRAY, JRCD, JMAC, JMATYP,
     &            MATLAYO, LACCFLA)
      
      IF (JRCD .EQ. 0) THEN
         PEEQ = ARRAY(1)
      ELSE
         PEEQ = 0.0D0
      END IF

      !-----------------------------------------------------------------
      ! 5. Assign Field Variable and Solution-Dependent State Variables
      !-----------------------------------------------------------------
      ! FIELD(1) controls interpolation in *PLASTIC, DEPENDENCIES=1
      FIELD(1) = DEPTH

      ! Store state variables for field contour visualization in ODB
      STATEV(1) = DEPTH
      STATEV(2) = PEEQ
      STATEV(3) = R_CURR

      RETURN
      END SUBROUTINE USDFLD
