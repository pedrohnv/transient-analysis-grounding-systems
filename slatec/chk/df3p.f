*DECK DF3P
      DOUBLE PRECISION FUNCTION DF3P (X)
C***BEGIN PROLOGUE  DF3P
C***PURPOSE  Subsidiary to
C***LIBRARY   SLATEC
C***AUTHOR  (UNKNOWN)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   ??????  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DF3P
      DOUBLE PRECISION X
C***FIRST EXECUTABLE STATEMENT  DF3P
      DF3P = 0.1D+01
      IF(X.GT.0.31415926535897932D+01) DF3P = 0.0D+00
      RETURN
      END