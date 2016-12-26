module ns

	implicit none

	include "fftw3.f"

	contains

	subroutine fourier_derivative(n, array)						! subroutine calculates derivative
											! in fourier space. 
		implicit none

		integer :: n, i
		double precision :: n1
		double complex, dimension(1:n) :: array

		do i = 1, n	 							! calculating derivative of input				
			if (i .le. n/2) then 						! function. This is done by multiplying
				n1 = i-1.0d0
				array(i) = array(i)*(0.0d0,1.0d0)*n1			! i*k (i is the imaginary unit) to every
			else								! element in the array. Here, k runs from
				n1 = (i-n-1)
				array(i) = array(i)*(0.0d0,1.0d0)*n1			! -l/2 to l/2 - 1. Note that the
			end if								! array is in wrapped around format only.
		end do
	
		array(n/2 + 1) = 0.0d0							! setting oddball wavenumber = 0

	end subroutine

	subroutine x_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(1:nx, 1:ny, 1:nz) :: input_3d
		double complex, dimension(1:nx) :: input_x
		integer*8 ::plan_x_fwd, plan_x_bwd
		integer :: fwd, bwd, deriv

		call dfftw_plan_dft_1d(plan_x_fwd, nx, input_x, input_x, fftw_forward, fftw_estimate)
		call dfftw_plan_dft_1d(plan_x_bwd, nx, input_x, input_x, fftw_backward, fftw_estimate)
		
		if (fwd .eq. 1) then
			do k = 1, nz
				do j = 1, ny										! transforming to fourier space.
					input_x = input_3d(:, j, k)							! copy to 1d array from 3d array.
					call dfftw_execute_dft(plan_x_fwd, input_x, input_x)				! call fftw 1d fft subroutine.
					if (deriv .eq. 1) then
						call fourier_derivative(nx, input_x)					! call derivative subroutine.
					end if
					input_3d(:, j, k) = input_x							! copy back to 3d input array from
				end do											! 1d array.
			end do
			input_3d = input_3d/nx										! renormalization after x-fft.
		end if

		if (bwd .eq. 1) then
			do k = 1, nz
				do j = 1, ny										! transforming back to physical space.
					input_x = input_3d(:, j, k)
					call dfftw_execute_dft(plan_x_bwd, input_x, input_x)
					input_3d(:, j, k) = input_x
				end do
			end do
		end if

		call dfftw_destroy_plan(plan_x_fwd)
		call dfftw_destroy_plan(plan_x_bwd)

	end subroutine

	subroutine y_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)

		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(1:nx, 1:ny, 1:nz) :: input_3d
		double complex, dimension(1:ny) :: input_y
		integer*8 ::plan_y_fwd, plan_y_bwd
		integer :: fwd, bwd, deriv

		call dfftw_plan_dft_1d(plan_y_fwd, ny, input_y, input_y, fftw_forward, fftw_estimate)
		call dfftw_plan_dft_1d(plan_y_bwd, ny, input_y, input_y, fftw_backward, fftw_estimate)
		
		if (fwd .eq. 1) then
			do k = 1, nz
				do i = 1, nx
					input_y = input_3d(i, :, k)
					call dfftw_execute_dft(plan_y_fwd, input_y, input_y)
					if (deriv .eq. 1) then
						call fourier_derivative(ny, input_y)
					end if
					input_3d(i, :, k) = input_y
				end do
			end do
			input_3d = input_3d/ny
		end if

		if (bwd .eq. 1) then
			do k = 1, nz
				do i = 1, nx
					input_y = input_3d(i, :, k)
					call dfftw_execute_dft(plan_y_bwd, input_y, input_y)
					input_3d(i, :, k) = input_y
				end do
			end do
		end if

		call dfftw_destroy_plan(plan_y_fwd)
		call dfftw_destroy_plan(plan_y_bwd)

	end subroutine

	subroutine z_derivative(input_3d, nx, ny, nz, fwd, bwd, deriv)
		
		implicit none

		integer :: i, j, k, nx, ny, nz		
		double complex, dimension(1:nx, 1:ny, 1:nz) :: input_3d
		double complex, dimension(1:nz) :: input_z
		integer*8 ::plan_z_fwd, plan_z_bwd
		integer :: fwd, bwd, deriv

		call dfftw_plan_dft_1d(plan_z_fwd, nz, input_z, input_z, fftw_forward, fftw_estimate)
		call dfftw_plan_dft_1d(plan_z_bwd, nz, input_z, input_z, fftw_backward, fftw_estimate)
		
		if (fwd .eq. 1) then 
			do j = 1, ny
				do i = 1, nx
					input_z = input_3d(i, j, :)
					call dfftw_execute_dft(plan_z_fwd, input_z, input_z)
					if (deriv .eq. 1) then
						call fourier_derivative(nz, input_z)
					end if
					input_3d(i, j, :) = input_z
				end do
			end do
			input_3d = input_3d/nz
		end if
		
		if (bwd .eq. 1) then
			do j = 1, ny
				do i = 1, nx
					input_z = input_3d(i, j, :)
					call dfftw_execute_dft(plan_z_bwd, input_z, input_z)
					input_3d(i, j, :) = input_z
				end do
			end do
		end if	

		call dfftw_destroy_plan(plan_z_fwd)
		call dfftw_destroy_plan(plan_z_bwd)

	end subroutine

	subroutine poisson(input, output, nx, ny, nz)

		integer :: nx, ny, nz, i, j, k
		double complex, dimension(1:nx,1:ny,1:nz) :: input, output

		call x_derivative(input, nx, ny, nz, 1, 0, 0)
		call y_derivative(input, nx, ny, nz, 1, 0, 0)
		call z_derivative(input, nx, ny, nz, 1, 0, 0)
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
					if ( i .eq. 1 .and. j .eq. 1 .and. k .eq. 1) then
						output(i,j,k) = 0.0
					else if (i .le. nx/2 .and. j .le. ny/2 .and. k .le. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-1)**2+(j-1)**2+(k-1)**2)
					else if (i .le. nx/2 .and. j .le. ny/2 .and. k .gt. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-1)**2+(j-1)**2+(k-nz-1)**2)
					else if (i .le. nx/2 .and. j .gt. ny/2 .and. k .le. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-1)**2+(j-ny-1)**2+(k-1)**2)
					else if (i .gt. nx/2 .and. j .le. ny/2 .and. k .le. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-nx-1)**2+(j-1)**2+(k-1)**2)
					else if (i .le. nx/2 .and. j .gt. ny/2 .and. k .gt. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-1)**2 + (j-ny-1)**2 + (k-nz-1)**2)
					else if (i .gt. nx/2 .and. j .gt. ny/2 .and. k .le. nz/2) then 
						output(i,j,k) = -input(i,j,k)/((i-nx-1)**2 + (j-ny-1)**2 + (k-1)**2)
					else if (i .gt. nx/2 .and. j .le. ny/2 .and. k .gt. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-nx-1)**2 + (j-1)**2 + (k-nz-1)**2) 
					else if (i .gt. nx/2 .and. j .gt. ny/2 .and. k .gt. nz/2) then
						output(i,j,k) = -input(i,j,k)/((i-nx-1)**2 + (j-ny-1)**2 + (k-nz-1)**2)
					end if
				end do
			end do
		end do
		call x_derivative(output, nx, ny, nz, 0, 1, 0)
		call y_derivative(output, nx, ny, nz, 0, 1, 0)
		call z_derivative(output, nx, ny, nz, 0, 1, 0)
		call x_derivative(input, nx, ny, nz, 0, 1, 0)
		call y_derivative(input, nx, ny, nz, 0, 1, 0)
		call z_derivative(input, nx, ny, nz, 0, 1, 0)		
					
	end subroutine

end module

program poissoncheck

	use ns

	implicit none

	integer, parameter :: nx=128, ny=128, nz=32
	integer :: i,j,k
	double complex, dimension(1:nx,1:ny,1:nz) :: rhs, result1
	double precision, dimension(1:nx) :: x
	double precision, dimension(1:ny) :: y
	double precision, dimension(1:nz) :: z
	double precision, parameter :: fx = 0, fy = 0, fz = 0, lx = 10, ly = 10, lz = 10, dx = (lx-fx)/nx, dy = (ly-fy)/ny, dz = (lz-fz)/nz
	double precision, parameter :: pi = acos(-1.0d0), a=1.0
	double precision, parameter :: x1 = 3.5, y1 = 5.0, x2 = 6.5, y2 = 5.0
!!$	double precision, parameter :: pi = acos(-1.0d0), dx = 2.0d0*pi/nx, dy = 2.0d0*pi/ny, dz = 2.0d0*pi/nz
!!$
!!$	do i = 1, nx
!!$		x(i) = (i-1)*dx
!!$	end do
!!$	do j = 1, ny
!!$		y(j) = (j-1)*dy
!!$	end do
!!$	do k = 1, nz
!!$		z(k) = (k-1)*dz
!!$	end do
!!$
!!$	do k = 1, nz
!!$		do j = 1, ny
!!$			do i = 1, nx	
!!$				rhs(i,j,k) = sin(5*x(i))*sin(5*y(j))*sin(8*z(k))
!!$				result(i,j,k) = 5*5*8*cos(5*x(i))*cos(5*y(j))*cos(8*z(k))
!!$			end do
!!$		end do
!!$	end do
!!$				
!!$	call x_derivative(rhs, nx, ny, nz, 1, 1, 1)
!!$	call y_derivative(rhs, nx, ny, nz, 1, 1, 1)
!!$	call z_derivative(rhs, nx, ny, nz, 1, 1, 1)
!!$
!!$	write(*,*) maxval(abs(real(rhs-result)))

	do i=1, nx
		x(i) = fx + (i-1)*dx
	end do
	do j=1,ny
		y(j) = fy + (j-1)*dy
	end do
	do k=1,nz
		z(k) = fz + (k-1)*dz
	end do

	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				! rhs(i,j,k) = exp(-10.0d0*(x(i)**2 + y(j)**2)) - (pi/(40.0*a**2))*(erf(sqrt(10.0)*a))**2 		Test Problem 1
				! rhs(i,j,k) = exp(-10.0d0*(2.0*x(i)**2 + 4.0*x(i)*y(j) + 5.0*y(j)**2)) - 0.128255/(4.0*a**2)		Test Problem 2
				  rhs(i,j,k) = exp(-((x(i)-3.5)**2 + (y(j)-5.0)**2 )/0.8) + exp(-((x(i)-6.5)**2 + (y(j)-5.0)**2)/0.8) - 5.02655/(100.0)
			end do
		end do
	end do

	write(*,*) maxval(real(rhs)), minval(real(rhs))

	call poisson(rhs, result1, nx, ny, nz)
	
	open(unit=25,file="result.dat")
        write(25,*) 'variables= "x","y","input","result"'
	write(25,*) 'zone I = 128 J = 128 F=POINT'
	do j = 1, ny
		do i = 1, nx
			write(25,'(4E30.18)') x(i), y(j), real(rhs(i,j,1)), real(result1(i,j,1))
		end do
	end do
	close(25)
end program
