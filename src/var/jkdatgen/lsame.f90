logical function lsame(ca, cb)

  implicit none

!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

!     .. Scalar Arguments ..
  character, intent(in) :: ca, cb

!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================

!     .. Local Scalars ..
  integer :: inta, intb, zcode

!     Test if the characters are equal

  lsame = ca .eq. cb
  if (lsame) return

!     Now test for equivalence if both characters are alphabetic.

  zcode = ichar('Z')

!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.

  inta = ichar(ca)
  intb = ichar(cb)

  if (zcode .eq. 90 .or. zcode .eq. 122) then

!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.

     if (inta .ge. 97 .and. inta .le. 122) inta = inta - 32
     if (intb .ge. 97 .and. intb .le. 122) intb = intb - 32

  else if (zcode .eq. 233 .or. zcode .eq. 169) then

!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.

     if (inta .ge. 129 .and. inta .le. 137 .or. &
          inta .ge. 145 .and. inta .le. 153 .or. &
          inta .ge. 162 .and. inta .le. 169) inta = inta + 64
     if (intb .ge. 129 .and. intb .le. 137 .or. &
          intb .ge. 145 .and. intb .le. 153 .or. &
          intb .ge. 162 .and. intb .le. 169) intb = intb + 64

  else if (zcode .eq. 218 .or. zcode .eq. 250) then

!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.

     if (inta .ge. 225 .and. inta .le. 250) inta = inta - 32
     if (intb .ge. 225 .and. intb .le. 250) intb = intb - 32
  end if
  lsame = inta .eq. intb

end function lsame
