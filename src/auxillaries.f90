MODULE AUXILLARIES
USE SECOND_PRECISION, ONLY : dp
USE CONSTANTS
IMPLICIT NONE
PUBLIC

CHARACTER(56) :: FMT_TIME    = '("+.",t18,82("."),t3,"Time: ",a,t100, "+")'
CHARACTER(56) :: FMT_CVU     = '("| ",t3, a,es10.3,t100, "|")'
CHARACTER(56) :: FMT10_CVU   = '("| ",t3, a,t13, es10.3,a,t100, "|")'
CHARACTER(56) :: FMT30_CVU   = '("| ",t3, a,t33, es10.3,a,t100, "|")'
CHARACTER(86) :: FMT10_2CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT10_3CVU  = '("| ",t3, a,t13, es10.3,a,t35,a,es10.3,a,t65,a,es10.3,a, t100, "|")'
CHARACTER(86) :: FMT_INTRMT  = '(" +",t4, 96("."), t4, a, t100, "+")'
CHARACTER(86) :: FMT_LEND    = '("+",t2, 98("."), t100, "+")'
CHARACTER(86) :: FMT_LOOPEND = '("+",t2, 49("--"), t100, "+")'
CHARACTER(86) :: FMT_SUB     = '("| ",t5,":",8("."), a,t100, "|")'
CHARACTER(86) :: FMT_MSG     = '("| ",t3,a,t100, "|")'
CHARACTER(86) :: FMT_HDR     = '("+",t2, 98(".") t3,a,t100, "+")'
CHARACTER(86) :: FMT_WARN0   = '("| WARNING:",t12,88("~"),"+", t12, a)'
CHARACTER(86) :: FMT_WARN1   = '("| WARNING:",t12,88("~"),"+", t12, a,f0.4)'
CHARACTER(86) :: FMT_NOTE0   = '("| NOTE:",t9,91("~"),"+", t9, a)'
CHARACTER(86) :: FMT_NOTE1   = '("| NOTE:",t9,91("~"),"+", t9, a,f0.4)'
CHARACTER(86) :: FMT_FAT0    = '("| FATAL ERROR:",t16,84("~"),"+", t16, a)'

CONTAINS

! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/cm3
real(dp) function C_AIR_cc(T,P)
  REAL(dp), INTENT(in) :: T, P
  C_AIR_cc = P/(kb*T)*1d-6
end function C_AIR_cc
! concentration of air. Pressure P in [Pa], temperature T in [K] returns molecules/m3
real(dp) function C_AIR_m3(T,P)
  REAL(dp), INTENT(in) :: T, P
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
    time = GTIME
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
  logical             :: invalue
  invalue = .false.
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
    READ(file_id, '(a)', iostat=ioi) buffer
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
! Function that linearily interpolates any value at current [time] (of
! type(timetype)) using timearray [conctime] (by default fractional day) and
! respective concentration timeseries [conc]. If [row] (integer) is provided,
! function will use that row as starting point for interpolation, or search for
! correct row if this fails. If optional unit is provided ('sec', 'min', 'hrs',
! 'day') then the corresponding time unit is used from [time].
!..............................................................................
REAL(dp) FUNCTION INTERP(conctime, conc, row, unit, timein)
  IMPLICIT NONE
  real(dp), OPTIONAL, intent(in)  :: timein
  REAL(dp), intent(in) :: conctime(:), conc(:)
  REAL(dp) :: x
  INTEGER, OPTIONAL :: row
  CHARACTER(*), OPTIONAL :: unit
  INTEGER :: rw, i
  if (PRESENT(timein)) THEN
    x = timein
  ELSE
    if (PRESENT(unit)) THEN
        if (unit .eq. 'sec') THEN
            x = GTIME%sec
        ELSEIF (unit .eq. 'min') THEN
            x = GTIME%min
        ELSEIF (unit .eq. 'hrs') THEN
            x = GTIME%hrs
        ELSEIF (unit .eq. 'day') THEN
            x = GTIME%day
        ELSE
            print FMT_WARN0, 'UNKNOWN TIME UNIT, can not interpolate, trying with days'
            x = GTIME%day
        END IF
    ELSE
        x = GTIME%day
    END IF
  END IF
  x = x + conctime(1)
  rw = 0
  if (PRESENT(row)) THEN
    rw = row
    if ((x > conctime(rw)) .and. (x < conctime(rw+1))) THEN
        continue
    else
      print FMT_WARN1,'Wrong row number is sent in to INTERP, searching for the real row x.', REAL(rw)
      rw = 0
      do WHILE (x>=conctime(rw+1))
        rw = rw + 1
      end do
      print FMT_NOTE1,'real row is: ', REAL(rw)
    end if
  else
    rw=1
    do i=1,size(conctime)
      if (conctime(i) > x) exit
    end do
    rw = i-1
    if (rw==size(conctime)) THEN
      INTERP = conc(rw)
      return
      ! rw = size(conctime)-1
    end if
  end if

    INTERP = (conc(rw+1)-conc(rw)) / (conctime(rw+1)-conctime(rw)) * (x-conctime(rw)) + conc(rw)
    if (NO_NEGATIVE_CONCENTRATIONS) INTERP = max(0d0, INTERP)

END FUNCTION INTERP

!==============================================================================
! Function finds out the internal index using the proper name. E.g. IndexFromName('APINENE')
! will return 150 in current version.
! Input:
! name: the name (trimmed or untrimmed) of the variable. This name must match the name in NAMES.DAT exactly.
! Output:
! integer
!..............................................................................
PURE INTEGER FUNCTION IndexFromName(NAME, list_of_names)
  IMPLICIT NONE
  character(*), INTENT(IN) :: NAME
  character(*), optional, INTENT(IN) :: list_of_names(:)
  integer :: i,m
  if (PRESENT(list_of_names)) then
    m = size(list_of_names)
    DO i=1, m
      if (TRIM(NAME) == TRIM(list_of_names(I))) EXIT
    END DO
  ELSE
    m = size(MODS)
    DO i=1, m
      if (TRIM(NAME) == TRIM(MODS(I)%NAME)) EXIT
    END DO
  END IF

  if (i == m+1) i=0
  IndexFromName = i

END FUNCTION IndexFromName

subroutine handle_file_io(ioi, file, halt)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ioi
  CHARACTER(*), INTENT(in) :: file
  CHARACTER(100) :: fmt
  CHARACTER(*), OPTIONAL, INTENT(in) :: halt
  if (PRESENT(halt)) THEN
    fmt = FMT_FAT0
  else
    fmt = FMT_WARN0
  end if
  if (ioi /= 0) THEN
    print fmt, 'Could not open '//TRIM(file)//', does it exist?'
    if (PRESENT(halt)) THEN
      print FMT_SUB, halt
      print FMT_LEND
      stop
    END IF
  END if
end subroutine handle_file_io

!====================================================================================
! f2chr and i2chr turns float and integer to characters, Scientific notation for real
!....................................................................................
PURE CHARACTER(LEN=12) FUNCTION f2chr(number)
    IMPLICIT NONE
    real(dp), INTENT(IN) :: number
    write(f2chr, '(es12.3)') number
    f2chr = ADJUSTL(f2chr)
END FUNCTION f2chr

PURE FUNCTION i2chr(number) result(out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: number
    CHARACTER(len=int(LOG10(MAX(ABS(number)*1d0, 1d0))+2)-min(0,sign(1,number))) :: out
    write(out, '(i0)') number
END FUNCTION i2chr

PURE FUNCTION di2chr(number) result(out)
    IMPLICIT NONE
    INTEGER(dint), INTENT(IN) :: number
    CHARACTER(len=int(LOG10(MAX(ABS(number)*1d0, 1d0))+2)-min(0,sign(1_dint,number))) :: out
    write(out, '(i0)') number
END FUNCTION di2chr

PURE LOGICAL FUNCTION equal(a,b)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a,b
    equal = ABS(a-b) .lt. 1d-200
END FUNCTION equal


!====================================================================================
! Calculates multimodal particle size distribution, used for intialization. Modevector
! is of length <number of modes>*3 and contains the ount median diameter, (CMD),
! geometric standard deviation (GSD) and relative size in total particle count
! (relative sizes must add up to 1). diameters is the diameter vector of the model
! (the bins), psd is the output vector and N is the total particle count in the output
! vector
!....................................................................................
subroutine Multimodal(modevector, diameters, psd, N)
    implicit none
    real(dp), INTENT(in)    :: modevector(:)
    real(dp), INTENT(in)    :: diameters(:)
    real(dp), INTENT(inout) :: psd(:)
    real(dp), allocatable   :: x(:), sumv(:)
    real(dp)                :: N ,mu,sig           ! size factor and total count
    INTEGER                 :: ii,nModes, uprInd, nb

    nb = size(psd)

    IF (MODULO(size(modevector),3) .ne. 0) THEN
        print*, 'The vector for modes is not correct'
        STOP
    ELSE
        nModes = size(modevector)/3
        if (GTIME%printnow) print*, 'Building PSD from ',nModes,' modes'
    END IF

    ALLOCATE(x(nb))
    ALLOCATE(sumv(nb))
    sumv = 0d0

    x = LOG10(diameters)

    DO ii = 1,nModes
        mu   = LOG10(modevector((ii-1)*3+1))
        sig  = modevector((ii-1)*3+2)
        sumv = sumv + modevector((ii-1)*3+3) * gauss(x, mu, sig)
    END DO

    psd = N * sumv / (sum( sumv ))
    uprInd = nb+1 - min(nb-1, 8)
    WHERE (psd(uprInd:) > 1d-12) psd(uprInd:) = 0d0

end subroutine Multimodal


pure function gauss(x,mu,sig) result(zz)
    real(dp), INTENT(IN) :: x(:)
    real(dp), INTENT(IN) :: mu,sig
    real(dp), allocatable :: zz(:)
    allocate(zz(size(x)))
    zz = exp(-(x-mu)**2/(2*sig**2))/sqrt(2*pi*sig**2)
end function gauss


pure elemental function saturation_conc_m3(A,B, Temperature) result(Vapour_concentration)
  real(dp), intent(in) :: A, B, temperature
  real(dp) :: Vapour_concentration, vapour_pressure

  ! Using antoine equation log_10(p) = A- (B/T)
  vapour_pressure      = 10 ** (A - (B/temperature)) ! in atm
  Vapour_concentration = (vapour_pressure*101325)/(kb * temperature) ! #/m3

end function saturation_conc_m3




end MODULE AUXILLARIES
