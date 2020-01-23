module acdc_system_extras

implicit none

contains

subroutine n_B_in_clusters(n_B)
	implicit none
	integer :: n_B(54)

	n_B = (/0,0,1,1,1,1,2,2,2,3,3,3,4,4,&
  5,5,0,0,0,0,0,1,1,1,1,1,2,2,2,&
  2,2,3,3,4,0,1,1,2,2,2,2,3,3,3,&
  3,3,4,4,4,5,1,1,0,0/)

end subroutine n_B_in_clusters
end module acdc_system_extras
