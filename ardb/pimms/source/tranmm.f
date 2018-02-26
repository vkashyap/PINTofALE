C+
	REAL FUNCTION TRANMM(E,COL)

        implicit none

        real e, col
C
C	RETURNS INTERSTELLAR TRANSMISSION USING BROWN & GOULD COEFFICIENTS
C	FOR PHOTON ENERGY 0.1 TO 10. KeV.
C	BELOW 0.1 KeV CROSS SECTIONS ARE EXTRAPOLATED (~E**(8/3))
C	ABOVE 10. KeV ABSORPTION < 1% FOR COLUMN DENSITIES <10**22 cm**2
C
C	ON INPUT:	E	= ENERGY (KEV)
C			COL	= COLUMN DENSITY (1E20/CM**2)
*                                                Now in plain /cm/cm
C
C	ON OUTPUT:	FUNCTION= INTERSTELLAR TRANSMISSION
C
C	N.B.	IN PROGRAM SPECFIT THIS ROUTINE IS CALLED VIA ENTRY POINT TRANIS
C		(THIS IS IN FILE PROFIL.MAC CREATED AT SPECFIT GENERATION TIME).
C-
C	AUTHOR:	PAUL LAMB (MSSL)			24-JUN-83
C       MODIFIED FOR MORRISON AND MCCAMMON XSECTIONS BY KOM 5-MAR-85
C	Modified:	 22-Jan-1985 Renamed procedure to TRANMM
C			 30-Apr-1986 Removed output of identification
C                        28-Sep-2001 Protected against underflow, KM
C
C       Modified 2003 September --- now calls TRANUVO if the photon
C                     energy is below 13.58eV (or lambda>912A)
C
	real TRANUVO

        real bgdif, e3, axs, arg
        integer i, ii
        real EX(14),C0(14),C1(14),C2(14)
        DATA EX/0.100,0.284,0.400,0.532,0.707,0.867,1.303,
     +          1.840,2.471,3.210,4.038,7.111,8.331,10.000/
        DATA C0/17.30,34.60,78.10,71.40,95.50,308.9,120.6,
     +          141.3,202.7,342.7,352.2,433.9,629.0,701.2/
        DATA C1/608.1,267.9,18.80,66.80,145.8,-380.6,169.3,
     +          146.8,104.7,18.70,18.70,-2.40,30.90,25.20/
        DATA C2/-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,
     +          -31.5,-17.0,0.000,0.000,0.750,0.000,0.000/
C
        IF(E.LT.0.)THEN
          E=0.
        ELSE IF(E.LT.0.01358)THEN
          TRANMM=TRANUVO(E,COL)
	ELSE
          DO 10 I=1,14
             II=I
          IF(E.LE.EX(I))GO TO 20
10        CONTINUE
20        E3=E*E*E
          AXS=(C0(II)+C1(II)*E+C2(II)*E*E)
        	BGDIF=AXS/E3
                BGDIF=-BGDIF*1.0E-4
          arg = col / 1.0e+20 * bgdif
*         ARG=COL*BGDIF
          TRANMM=0.
          IF(ARG.GT.-60.0)TRANMM=EXP(ARG)
        END IF
        END
