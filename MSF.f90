program MSF
implicit none
	
	real(8) :: xx, yy, zz ,kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kz1, kz2, kz3, kz4, h
	real(8), dimension(3,3) :: ko1, ko2, ko3, ko4, o
	real(8), dimension(3,3) :: DF	!Jacobian Matrix of the system function
	real(8), dimension(3,3) :: DH
	integer :: i, j, ijk, j0, j1, j2, it	!counter
	real(8), external :: f	!function to calculate
	real(8), parameter :: a = 10., b = 28., c = 2.	!parameter of the function
	real(8), dimension(3) :: lyn, cc, dd
	real(8), dimension(3,3) :: Ot, QQ, RR, Rt, Qt
	real(8) :: temp_Lya, p
	integer ::  time0, time1, q
	LOGICAL sing
	
	call system_clock(time0)
	
	open (11, file = '1_1.dat')
	open (12, file = '1_2.dat')
	open (13, file = '1_3.dat')
	open (21, file = '2_1.dat')
	open (22, file = '2_2.dat')
	open (23, file = '2_3.dat')
	open (31, file = '3_1.dat')
	open (32, file = '3_2.dat')
	open (33, file = '3_3.dat')
 
	h=0.001d0
	
	do q = 1, 9
		DH = 0	
		do p = 0, 300, 2
		
			select case (q)
				case(1)
					DH(1, 1) = p / 3.0d0 
				case(2)
					DH(1, 2) = p / 3.0d0
				case(3)
					DH(1, 3) = p / 3.0d0
				case(4)
					DH(2, 1) = p / 3.0d0
				case(5)
					DH(2, 2) = p / 3.0d0
				case(6)
					DH(2, 3) = p / 3.0d0
				case(7)
					DH(3, 1) = p / 3.0d0
				case(8)
					DH(3, 2) = p / 3.0d0
				case(9)
					DH(3, 3) = p / 3.0d0
			end select

!			Jacobian Matrix	
			DF(1, 1) = -a - DH(1, 1)
			DF(1, 2) = a - DH(1, 2)
			DF(1, 3) = 0 - DH(1, 3)
			DF(2, 2) = -1 - DH(2, 2)
			DF(3, 3) = -c - DH(3, 3)
			
			ijk = 0
!			define initial value
			call random_seed()
			call random_number(xx)
			call random_number(yy)
			call random_number(zz)
			o = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3, 3/))
			QQ = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3, 3/))	
			RR = 0
			Ot = 0
			Rt = 0
			Qt = 0
			
			
		!	Solve equations using RK4	
			do i = 1, 5000000

				kx1 = f(1, xx, yy, zz, a, b, c)
				ky1 = f(2, xx, yy, zz, a, b, c)
				kz1 = f(3, xx, yy, zz, a, b, c)	
				if (i > 1000000) then		
					DF(2, 1) = b - zz - DH(2, 1)
					DF(2, 3) = -xx - DH(2, 3)
					DF(3, 1) = yy - DH(3, 1)
					DF(3, 2) = xx - DH(3, 2)
					ko1 = matmul(DF, o)
				end if
				
				kx2 = f(1, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
				ky2 = f(2, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
				kz2 = f(3, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
				if (i > 1000000) then	
					DF(2, 1) = b - (zz + kz1 * h / 2) - DH(2, 1)
					DF(2, 3) = -(xx + kx1 * h / 2) - DH(2, 3)
					DF(3, 1) = yy + ky1 * h / 2 - DH(3, 1)
					DF(3, 2) = xx + kx1 * h / 2 - DH(3, 2)
					ko2 = matmul(DF, o + ko1 * h / 2)
				end if 
				
				kx3 = f(1, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
				ky3 = f(2, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
				kz3 = f(3, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
				if (i > 1000000) then		
					DF(2, 1) = b - (zz + kz2 * h / 2) - DH(2, 1)
					DF(2, 3) = -(xx + kx2 * h / 2) - DH(2, 3)
					DF(3, 1) = yy + ky2 * h / 2 - DH(3, 1)
					DF(3, 2) = xx + kx2 * h / 2 - DH(3, 2)
					ko3 = matmul(DF, o + ko2 * h / 2)
				end if
				
				kx4 = f(1, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
				ky4 = f(2, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
				kz4 = f(3, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
				if (i > 1000000) then	
					DF(2, 1) = b - (zz + kz3 * h) - DH(2, 1)
					DF(2, 3) = -(xx + kx3 * h) - DH(2, 3)
					DF(3, 1) = yy + ky3 * h - DH(3, 1)
					DF(3, 2) = xx + kx3 * h - DH(3, 2)
					ko4 = matmul(DF, o + ko3 * h)
				end if

				xx = xx + h * (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6
				yy = yy + h * (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6
				zz = zz + h * (kz1 + 2 * kz2 + 2 * kz3 + kz4) / 6

				if (i > 1000000) then
					ijk = ijk + 1
					o = o + h * (ko1 + 2 * ko2 + 2 * ko3 + ko4) / 6

					Ot = matmul(o, QQ)
					
		!			QR decompostion
					call ddqrdcmp(Ot, 3, cc, dd, sing)
					
					do j=1, 2
						Rt(j,j) = dd(j)
						do j1 = j+1, 3
							Rt(j, j1) = Ot(j, j1)
							Ot(j, j1) = 0
						end do
					end do
					j = 3
					Rt(j, j) = dd(j)
					
					QQ = 0
					do j=1,3
						QQ(j,j)=1
					end do
					
					do j = 1, 2
						do j1=1, 3
							do j2=1, 3
								Qt(j1, j2) = -Ot(j1, j) * Ot(j2, j) / cc(j)
								if (j1 == j2) then
									Qt(j1, j2) = Qt(j1, j2)+1
								end if
							end do
						end do
						QQ = matmul(QQ, Qt)
					end do

					do j=1, 3
						RR(j, j)=RR(j, j)+log(abs(Rt(j, j)))
					end do

		!			reset o
					o = 0
					do it= 1, 3
						o(it, it) = 1
					end do

				end if
				
			end do
			
			do j = 1, 3
				lyn(j) = RR(j, j)
			end do
			temp_Lya = max(lyn(1), lyn(2), lyn(3))
			temp_Lya = temp_Lya / ijk / h
		
			select case (q)
				case(1)
					write(11, *), DH(1, 1), temp_Lya 
				case(2)
					write(21, *), DH(1, 2), temp_Lya 
				case(3)
					write(31, *), DH(1, 3), temp_Lya 
				case(4)
					write(12, *), DH(2, 1), temp_Lya 
				case(5)
					write(22, *), DH(2, 2), temp_Lya 
				case(6)
					write(32, *), DH(2, 3), temp_Lya 
				case(7)
					write(13, *), DH(3, 1), temp_Lya 
				case(8)
					write(23, *), DH(3, 2), temp_Lya 
				case(9)
					write(33, *), DH(3, 3), temp_Lya 
			end select			
		end do
		
		call system_clock(time1)
		write(*,*) 'Elapsed time: ', (time1 - time0)*1E-4
		
	end do
	
end program MSF





	function f(j, x, y, z, a, b, c)
	implicit none
		real(8) :: f
		real(8) :: x, y, z
		real(8) :: a, b, c
		integer :: j
		if (j == 1) then
			f = a * (y - x)
		else if (j == 2) then
			f = x * (b - z) - y
		else if (j == 3) then
			f =x * y - c * z
		end if
	end function f




	SUBROUTINE ddqrdcmp(a, n, c, d, sing)
	implicit none
		INTEGER n,np
		REAL(8) :: a(n, n), c(n), d(n)
		LOGICAL :: sing
		INTEGER :: i, j, k
		REAL(8) :: scale, sigma, sum, tau

		sing = .false.

		do k = 1, n-1
			scale = 0.
			do i = k, n
				scale = max(scale, abs(a(i, k)))
			end do

			if(scale == 0.)then
				sing = .true.
				c(k) = 0.
				d(k) = 0.
			else
				do i = k, n
					a(i, k) = a(i, k) / scale
				end do
				sum = 0.
				do i=k, n
					sum = sum + a(i,k) ** 2
				end do
				sigma = sign(sqrt(sum), a(k, k))
				a(k, k) = a(k, k) + sigma
				c(k) = sigma * a(k, k)
				d(k) = -scale * sigma
				do j=k + 1, n
					sum = 0.
					do i = k, n
						sum = sum + a(i, k) * a(i, j)
					end do
					tau = sum / c(k)
					do i = k, n
						a(i, j) = a(i, j) - tau * a(i, k)
					end do
				end do
			end if
		end do

		d(n)=a(n, n)

		if(d(n) == 0.) then
			sing=.true.
		end if

	END