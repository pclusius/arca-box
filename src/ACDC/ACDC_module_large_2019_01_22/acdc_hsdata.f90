module acdc_hsdata

implicit none

contains

subroutine get_hsdata(molec_in,lhs,deltah_out,deltas_out)
	implicit none
	integer, intent(in) :: molec_in(2)						! Molecule numbers in the order A, N
	logical, intent(out) :: lhs								! Logical telling if data exists for the given composition
	real(kind(1.d0)), intent(out) :: deltah_out, deltas_out ! Thermochemistry (kcal/mol, cal/mol/K)
	
	type hsdata_type
		integer:: molec(2)
		real(kind(1.d0)) :: deltah, deltas
	end type hsdata_type
	type(hsdata_type), save :: qcdata(16)
	integer, save :: ntot = 16
	integer :: n
	logical, save :: firstcall = .true.
	
	if (firstcall) then
		firstcall = .false.
		n = 0
		
		! RICC2//B3LYP data
		n = n+1
		qcdata(n)%molec(1:2) = (/1, 0/)
		qcdata(n)%deltah = 0.d0
		qcdata(n)%deltas = 0.d0
		
		n = n+1
		qcdata(n)%molec(1:2) = (/0, 1/)
		qcdata(n)%deltah = 0.d0
		qcdata(n)%deltas = 0.d0
		
		n = n+1
		qcdata(n)%molec(1:2) = (/2, 0/)
		qcdata(n)%deltah = -17.8487481487d0
		qcdata(n)%deltas = -33.418352d0

		n = n+1
		qcdata(n)%molec(1:2) = (/1, 1/)
		qcdata(n)%deltah = -15.99849325487d0
		qcdata(n)%deltas = -28.136935d0

		n = n+1
		qcdata(n)%molec(1:2) = (/2, 1/)
		qcdata(n)%deltah = -44.9962620190d0
		qcdata(n)%deltas = -71.024601d0

		n = n+1
		qcdata(n)%molec(1:2) = (/3, 1/)
		qcdata(n)%deltah = -66.0568351924d0
		qcdata(n)%deltas = -107.724444d0

		n = n+1
		qcdata(n)%molec(1:2) = (/2, 2/)
		qcdata(n)%deltah = -64.4605106074d0
		qcdata(n)%deltas = -104.450256d0

		n = n+1
		qcdata(n)%molec(1:2) = (/3, 2/)
		qcdata(n)%deltah = -92.0891545072d0
		qcdata(n)%deltas = -143.176618d0

		n = n+1
		qcdata(n)%molec(1:2) = (/4, 2/)
		qcdata(n)%deltah = -115.1292995691d0
		qcdata(n)%deltas = -183.339004d0

		n = n+1
		qcdata(n)%molec(1:2) = (/3, 3/)
		qcdata(n)%deltah = -116.59903394959d0
		qcdata(n)%deltas = -177.985175d0

		n = n+1
		qcdata(n)%molec(1:2) = (/4, 3/)
		qcdata(n)%deltah = -145.1669199422d0
		qcdata(n)%deltas = -222.332473d0

		n = n+1
		qcdata(n)%molec(1:2) = (/5, 3/)
		qcdata(n)%deltah = -168.7911001831d0
		qcdata(n)%deltas = -260.547804d0

		n = n+1
		qcdata(n)%molec(1:2) = (/4, 4/)
		qcdata(n)%deltah = -164.3486235509d0
		qcdata(n)%deltas = -251.025092d0

		n = n+1
		qcdata(n)%molec(1:2) = (/5, 4/)
		qcdata(n)%deltah = -191.8609959824d0
		qcdata(n)%deltas = -291.053716d0

		n = n+1
		qcdata(n)%molec(1:2) = (/4, 5/)
		qcdata(n)%deltah = -186.4720560274d0
		qcdata(n)%deltas = -296.507159d0

		n = n+1
		qcdata(n)%molec(1:2) = (/5, 5/)
		qcdata(n)%deltah = -221.6506206689d0
		qcdata(n)%deltas = -332.489463d0

	end if
	
	lhs = .false.
	deltah_out = 0.d0
	deltas_out = 0.d0
	do n = 1,ntot
		if (all(molec_in .eq. qcdata(n)%molec)) then
			lhs = .true.
			deltah_out = qcdata(n)%deltah
			deltas_out = qcdata(n)%deltas
			return
		end if
	end do
	
end subroutine get_hsdata


end module acdc_hsdata