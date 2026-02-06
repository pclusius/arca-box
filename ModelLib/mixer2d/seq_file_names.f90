module seq_file_names
! program seq_file_names
implicit None

contains

integer function findex(prefix,n,suffix) result(i)
implicit None

character(len=*), intent(in)  :: prefix,suffix
integer, intent(in)           :: n
character(len=128)            :: fname
character(len=10)             :: fmt
write(fmt,'(a,i0,a)') '(a,i0.',n,',a)'

i=0
do
  write(fname, fmt) prefix,i,suffix
  if (access(TRIM(fname),' ') /= 0) exit
  i = i + 1
end do
i = i-1

end function findex

end module seq_file_names
! end program seq_file_names
