MODULE AUXILLARIES
USE SECOND_PRECISION, ONLY : dp
USE CONSTANTS
IMPLICIT NONE
PUBLIC

CHARACTER(56) :: FMT_TIME    = '("+.",t18,82("."),t3,"Time: ",a,t100, "+")'
CHARACTER(56) :: FMT10_CVU   = '("| ",t3, a,t13, es10.3,a,t100, "|")'
CHARACTER(56) :: FMT30_CVU   = '("| ",t3, a,t33, es10.3,a,t100, "|")'
CHARACTER(86) :: FMT10_2CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT10_3CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a,t65,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT_LEND    = '("+",t2, 98(".") t100, "+")'
CHARACTER(86) :: FMT_SUB     = '("| ",t5,":",8(".") a,t100, "|")'
CHARACTER(86) :: FMT_MSG     = '("| ",t3,a,t100, "|")'
CHARACTER(86) :: FMT_HDR     = '("+",t2, 98(".") t3,a,t100, "+")'
CHARACTER(86) :: FMT_WARN0   = '("| WARNING:",t12,88("~"),"+", t12, a)'
CHARACTER(86) :: FMT_WARN1   = '("| WARNING:",t12,88("~"),"+", t12, a,f0.4)'
CHARACTER(86) :: FMT_NOTE0   = '("| NOTE:",t9,91("~"),"+", t9, a)'
CHARACTER(86) :: FMT_NOTE1   = '("| NOTE:",t9,91("~"),"+", t9, a,f0.4)'
CHARACTER(86) :: FMT_FAT0   = '("| FATAL ERROR:",t16,84("~"),"+", t16, a)'

CONTAINS

! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/cm3
real(dp) function C_AIR_cc(T,P)
  REAL(dp) :: T, P
  C_AIR_cc = P/(kb*T)*1d-6
end function C_AIR_cc
! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/m3
real(dp) function C_AIR_m3(T,P)
  REAL(dp) :: T, P
  C_AIR_m3 = P/(kb*T)
end function C_AIR_m3

! Converts hours to seconds
real(dp) function hrs_to_s(h)
  REAL(dp) :: h
  hrs_to_s = h*3600
end function hrs_to_s

real(dp) function sec_to_d(s)
  REAL(dp) :: s
  sec_to_d = s/24d0/3600d0
end function sec_to_d

! Checks if model time is at some exact minute interval
! t = time in seconds (real8)
! check = minute that is to be checked (integer)
logical function EVENMIN(t,check,zero)
  IMPLICIT NONE
  INTEGER :: check, i
  INTEGER, optional :: zero
  REAL(dp):: t
  i = 0
  if (PRESENT(zero)) i=-1
  IF ((MODULO(int(t+0.5d0), 60*check) == 0) .and. (int(t+0.5d0)>i)) THEN
    EVENMIN = .true.
  ELSE
    EVENMIN = .false.
  END IF
end function EVENMIN

!==============================================================================
! Function to return y-value at [time] from a sine function. Use to model a
! periodic variable such as temperature or relative humidity. By default uses
! 24 hour period.
!..............................................................................
REAL(dp) FUNCTION PERIODICAL(MODS, timein)
  IMPLICIT NONE
  type(input_mod)   :: MODS
  type(timetype), OPTIONAL  :: timein
  type(timetype)            :: time
  real(dp)                  :: A
  if (PRESENT(timein)) THEN
    time = timein
  ELSE
    time = MODELTIME
  END IF
  A = (MODS%max-MODS%min)/2d0
  PERIODICAL = A*SIN(MODS%fv*2d0*pi*time%sec/time%SIM_TIME_S &
              - 2d0*pi*MODS%ph/time%SIM_TIME_H) + MODS%min + A
end FUNCTION PERIODICAL

!==============================================================================
! Function to count number of rows in a file. Takes in the UNIT (file pointer)
! and returns number of rows in a file (integer). The function can be used at
! any point of file reading, the UNIT will be moved back where it was in the
! beginning.
! NOTE:
! The function does not check if the file has equal amounts of colums in each row
! Empty rows are also counted
!..............................................................................
INTEGER FUNCTION ROWCOUNT(file_id,hdr)
  IMPLICIT NONE
  INTEGER :: file_id, ioi
  INTEGER(8) :: OFFSET
  character(1), optional :: hdr
  character(20) :: buffer
  OFFSET = FTELL(file_id)
  REWIND(file_id)
  ioi = 0
  ROWCOUNT = 0
  DO WHILE (ioi==0)
    READ(file_id, *, iostat=ioi) buffer
      if (present(hdr)) THEN
        if (trim(buffer(1:1)) .ne. hdr) ROWCOUNT = ROWCOUNT + 1
      ELSE
        ROWCOUNT = ROWCOUNT + 1
      END IF
  END DO
  ROWCOUNT = ROWCOUNT-1
  REWIND(file_id)
  CALL FSEEK(file_id, OFFSET, 0)
END FUNCTION ROWCOUNT


!==============================================================================
! Function to count number of colums in a file. Takes in the
! UNIT (file pointer) and returns number of columns in the FIRST row in file
! (integer). The function can be used at any point of file reading, the UNIT
! will be moved back where it was in the beginning.
! NOTE:
! Optionally separator can be sent in, e.g. for tab: COLCOUNT(file_id, achar(9))
! This function will still count all consecutive tabs as 1!
! The function does not check if the file has equal amounts of colums in each row
! Empty row or row with only spaces will return 0
!..............................................................................
INTEGER FUNCTION COLCOUNT(file_id, separator)
  IMPLICIT NONE
  INTEGER             :: file_id,ioi,i
  INTEGER(8)          :: OFFSET
  CHARACTER(6000)     :: buffer
  CHARACTER,OPTIONAL  :: separator
  CHARACTER           :: sep
  logical             :: invalue = .false.

  if(present(separator)) then
    sep = separator
  else
    sep = " "
  end if
  ioi = 0
  COLCOUNT = 0
  OFFSET = FTELL(file_id)
  REWIND(file_id)
  do i=1,2
    READ(file_id, '(a)', iostat=ioi), buffer
  end do
  if (ioi /= 0) THEN
    COLCOUNT = -9999
  else
    DO i=1,LEN_TRIM(buffer)
      if ((buffer(i:i) /= sep) .and. (invalue .neqv. .true.)) THEN
        if (COLCOUNT == 0) COLCOUNT = 1
        invalue = .true.
      elseif ((buffer(i:i) == sep) .and. (invalue) ) THEN
        invalue = .false.
        COLCOUNT = COLCOUNT +1
      end if
    END DO
  end if

  REWIND(file_id)
  CALL FSEEK(file_id, OFFSET, 0)
END FUNCTION COLCOUNT

!==============================================================================
! Gets the current row in timearray [conctime] for which:
! [conctime(getrow)] < time < [conctime(getrow+1)]
! Can be used if for some reason we have different time intervals for different
! inputs. The same  functionality is also incorporated in INTERP()
!..............................................................................
REAL(dp) FUNCTION getrow(time, conctime)
  IMPLICIT NONE
  INTEGER :: row = 0
  REAL(dp) :: time, conctime(:)
  conctime = conctime - conctime(1)
  do WHILE (time/3600d0/24d0>=conctime(row+1))
    row = row + 1
  end do
  getrow = row
END FUNCTION getrow

!==============================================================================
! Function that linearily interpolates any value at current [time] (of
! type(timetype)) using timearray [conctime] (by default fractional day) and
! respective concentration timeseries [conc]. If [row] (integer) is provided,
! function will use that row as starting point for interpolation, or search for
! correct row if this fails. If optional unit is provided ('sec', 'min', 'hrs',
! 'day') then the corresponding time unit is used from [time].
!..............................................................................
REAL(dp) FUNCTION INTERP(conctime, conc, row, unit, timein)
  IMPLICIT NONE
  type(timetype), OPTIONAL  :: timein
  type(timetype)  :: time
  REAL(dp) :: conctime(:), conc(:), now
  INTEGER, OPTIONAL :: row
  CHARACTER(*), OPTIONAL :: unit
  INTEGER :: rw
  if (PRESENT(timein)) THEN
    time = timein
  ELSE
    time = MODELTIME
  END IF
  if (PRESENT(unit)) THEN
    if (unit .eq. 'sec') THEN
      now = time%sec
    ELSEIF (unit .eq. 'min') THEN
      now = time%min
    ELSEIF (unit .eq. 'hrs') THEN
      now = time%hrs
    ELSEIF (unit .eq. 'day') THEN
      now = time%day
    ELSE
      print FMT_WARN0, 'UNKNOWN TIME UNIT, can not interpolate, trying with days'
      now = time%day
    END IF
  ELSE
    now = time%day
  END IF

  conctime = conctime - conctime(1)
  rw = 0
  if (PRESENT(row)) THEN
    rw = row
    if ((now > conctime(rw)) .and. (now < conctime(rw+1))) THEN
        continue
    else
      print FMT_WARN1,'Wrong row number is sent in to INTERP, searching for the real row now.', REAL(rw)
      rw = 0
      do WHILE (now>=conctime(rw+1))
        rw = rw + 1
      end do
      print FMT_NOTE1,'real row is: ', REAL(rw)
    end if
  else
    do WHILE (now>=conctime(rw+1))
      rw = rw + 1
    end do

  end if
  INTERP = (now-conctime(rw))/(conctime(rw+1)-conctime(rw))*(conc(rw+1)-conc(rw))+conc(rw)
END  FUNCTION INTERP


INTEGER FUNCTION CNTNONTYPESNML(NMLSTR)
  IMPLICIT NONE
  character(*):: NMLSTR
  INTEGER I, types
  types = 0
  CNTNONTYPESNML = 0
  DO I=1,LEN(TRIM(NMLSTR))
    IF (NMLSTR(I:I) == '=') CNTNONTYPESNML = CNTNONTYPESNML + 1
    IF (NMLSTR(I:I) == '%') types = types + 1
  END DO
  CNTNONTYPESNML = CNTNONTYPESNML - types
END FUNCTION CNTNONTYPESNML

INTEGER FUNCTION IndexFromName(name)
  IMPLICIT NONE
  character(*) :: name
  integer :: i
  DO i=1, size(MODS)
    if (TRIM(NAME) == TRIM(MODS(I)%NAME)) EXIT
  END DO
  IndexFromName = i
END FUNCTION IndexFromName

end MODULE AUXILLARIES
