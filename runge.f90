!> \brief An embedded Runge-Kutta solver

!> This subroutine solves ordinary differential equations
!> using the embedded RK routine whose coefficients are
!> stored in \c params .
!> \param n The number of dependent variables
!> \param nstep Maximum number of iterations to take
!> \param tol Error tolerance
!> \param tinit Initial value of the independent parameter
!> \param hinit Starting integration step size
!> \param hmax Maximum step size
!> \param upf Multiplicative parameter deciding percentage of current \c h value which becomes the new \c h'
!> \param x Location in state space
!> \param grad Pointer to function which returns the righthand or forcing function (derivative)
!> \param term Pointer to function which decides whether to terminate on each step
!> \param useterm Integer either 1 or 0 indicating either use or do not use \c term to decide to terminate
!> \param out Pointer to a function which outputs data at each step
!> \param useout Integer either 1 or 0 indicating either to use or do not use \c out to output data
!> \param data Customizable data to pass to the \c grad subroutine
!> \param params File containing coefficients defining the embedded Runge-Kutta method to be used
!> \param vbosel Verbose level; 0 for no output, > 0 for successively more output
subroutine dgerkem(n,nstep,tol,tinit,hinit,hmax,upf,x,grad,term,useterm,out,useout,data,params,vbosel)
	implicit none
	integer :: n
	integer :: nstep
	real*8 :: tol,tinit,hinit,hmax,upf
	real*8 :: x(n)
	external grad
	external term
	integer :: useterm
	external out
	integer :: useout
	external data
	character (256) :: params
	integer :: vbosel

	! Local variables
	integer :: i,j,k,m,res
	integer :: nstage,order1,order2
	real*8 :: t,tt,h,emin,s,tot
	real*8, allocatable :: acoef(:),ccoef(:),dcoef(:)
	real*8, allocatable :: bcoef(:,:),kv(:,:)
	real*8, allocatable :: cur(:),drv(:),ev(:),xt(:)

	! Load parameters from files
	open (1,FILE=params)
	read (1,*) nstage,order1,order2
	allocate(acoef(nstage),bcoef(nstage,nstage),ccoef(nstage),dcoef(nstage),kv(n,nstage),cur(n),drv(n),ev(n),xt(n))
	read (1,*) acoef ! The timestep coefficients
	read (1,*) bcoef ! The stage-wise intermediate coefficients
	read (1,*) ccoef ! Lower order set of evaluation coefficients
	read (1,*) dcoef ! Higher order set of evaluation coefficients, embedded by definition
	close(1)

	! Start verbosity
	if ( vbosel > 0 ) then
		write (*,*) "----------------------------------------------------------------------------------------------------------"
		write (*,"(A8)",advance="no") 'i'
		write (*,"(A16)",advance="no") 't'
		write (*,"(A16)",advance="no") 'h'
		write (*,"(A16)",advance="no") 'err'
		write (*,"(A16)",advance="no") "h'"
		write (*,"(A16)",advance="no") "diff"
		write (*,"(A16)") 'hmax'
		write (*,*) "----------------------------------------------------------------------------------------------------------"
	endif

	! Start the process
	t = tinit
	h = hinit

	! Write out initial position
	if ( vbosel > 0 ) then
		write (*,"(8I8)",advance="no") 0
		write (*,"(F16.6)",advance="no") t
		write (*,"(F16.6) ",advance="no") h
		write (*,"(A16)",advance="no") ''
		write (*,"(A16)",advance="no") ''
		write (*,"(A16)",advance="no") ''
		write (*,"(A16)") ''
	endif

	! Start outer iteration
	do i = 1,nstep
		! Clear the k-vector kv
		do j = 1,n
			do k = 1,nstage
				kv(j,k) = 0.0
			enddo
		enddo
		do j = 1,nstage
			do k = 1,n
				cur(k) = x(k)
			enddo
			do k = 1,j
				do m = 1,n
					cur(m) = cur(m) + bcoef(k,j) * kv(m,k)
				enddo
			enddo
			tt = t + acoef(j) * h
			call grad (n,tt,cur,drv,data)
			do k = 1,n
				kv(k,j) = h * drv(k)
			enddo
		enddo

		! Calculate the error using embedded formula
		do j = 1,n
			ev(j) = 0.0
		enddo
		do j = 1,n
			do k = 1,nstage
				ev(j) = ev(j) + ( dcoef(k) - ccoef(k) ) * kv(j,k)
			enddo
			ev(j) = abs( tol / ev(j) ) ** ( 1.0 / real( order2 ) )
		enddo

		! Calculate the minimum error in the vector
		do j = 1,n
			if ( ev(j) < emin .OR. j == 1 ) then
				emin = ev(j)
			endif
		enddo

		! Use emin to calculate the new stepsize h
		s = upf * emin * h
		if ( s < hmax ) then
			h = s
		else
			h = hmax
		endif

		! Clear the k-vector kv, again
		do j = 1,n
			do k = 1,nstage
				kv(j,k) = 0.0
			enddo
		enddo

		! Perform the Runge-Kutta calculation again using new h
		do j = 1,nstage
			do k = 1,n
				cur(k) = x(k)
			enddo
			do k = 1,j
				do m = 1,n
					cur(m) = cur(m) + bcoef(k,j) * kv(m,k)
				enddo
			enddo
			tt = t + acoef(j) * h
			call grad (n,tt,cur,drv,data)
			do k = 1,n
				kv(k,j) = h * drv(k)
			enddo
		enddo
		do j = 1,n
			xt(j) = x(j) ! Copy then modify
			do k = 1,nstage
				x(j) = x(j) + ccoef(k) * kv(j,k);
			enddo
		enddo
		tt = t
		t = t + h

		if ( vbosel > 0 ) then
			tot = 0.0
			do j = 1,n
				tot = tot + ( xt(j) - x(j) )**2
			enddo
			tot = sqrt (tot)
			write (*,"(8I8)",advance="no") i
			write (*,"(F16.6)",advance="no") t
			write (*,"(F16.6) ",advance="no") h
			write (*,"(F16.6)",advance="no") emin
			write (*,"(F16.6)",advance="no") s
			write (*,"(F16.6)",advance="no") tot
			write (*,"(F16.6)") hmax
		endif

		if ( useout .NE. 0 ) then
			call out (n,t,tt,x,xt,drv,data,res)
			if ( res .NE. 0 ) then
				! Error routine
			endif
		endif

		if ( useterm .NE. 0 ) then
			call term (n,t,tt,x,xt,drv,data,res)
			if ( res .NE. 0 ) then ! res > 0 implies termination
				exit
			endif
		endif
	enddo
	deallocate(acoef,bcoef,ccoef,dcoef,kv,cur,drv,ev,xt)
end subroutine

!> Version of \c dgerkem which takes RK parameters in arrays directly and takes one step per call

!> This is the same function as \c dgerkem with the exception that
!> instead of taking the file \c params containing the RK parameters,
!> it takes the values \c nstage , \c order , \c acoef , \c bcoef ,
!> \c ccoef , \c dcoef , and it only performs a single step in the
!> RK scheme.
!> \param n The number of dependent variables
!> \param t The independent parameter
!> \param h Starting integration step size
!> \param hmax Maximum step size
!> \param tol Error tolerance
!> \param upf Multiplicative parameter deciding percentage of current \c h value which becomes the new \c h'
!> \param x Location in state space - input
!> \param xt Location in state space - output
!> \param grad Pointer to function which returns the righthand or forcing function (derivative)
!> \param data Customizable data to pass to the \c grad subroutine
!> \param nstage Number of stages in the embedded Runge-Kutta scheme
!> \param order1 The roundoff error for the lower order scheme
!> \param order2 The roundoff error for the higher order scheme
!> \param acoef The time step parameters
!> \param bcoef Determines the values at which the dependent parameters are evaluated
!> \param ccoef Lower order multipliers on the \c kv values
!> \param dcoef Higher order multiplier on the \c kv values
subroutine dgerkems(n,t,h,hmax,tol,upf,x,xt,grad,data,nstage,order1,order2,acoef,bcoef,ccoef,dcoef)
	implicit none
	integer :: n
	real*8 :: t,h,hmax,tol,upf
	real*8 :: x(n),xt(n)
	external grad
	external data
	integer :: nstage,order1,order2
	real*8 :: acoef(nstage),bcoef(nstage,nstage),ccoef(nstage),dcoef(nstage)

	! Local variables
	integer :: j,k,m,res
	real*8 :: tt,emin,s
	real*8 :: kv(n,nstage)
	real*8 :: drv(n),ev(n)

	! Clear the k-vector kv
	do j = 1,n
		do k = 1,nstage
			kv(j,k) = 0.0
		enddo
	enddo
	do j = 1,nstage
		do k = 1,n
			xt(k) = x(k)
		enddo
		do k = 1,j
			do m = 1,n
				xt(m) = xt(m) + bcoef(k,j) * kv(m,k)
			enddo
		enddo
		tt = t + acoef(j) * h
		call grad (n,tt,xt,drv,data)
		do k = 1,n
			kv(k,j) = h * drv(k)
		enddo
	enddo

	! Calculate the error using embedded formula
	do j = 1,n
		ev(j) = 0.0
	enddo
	do j = 1,n
		do k = 1,nstage
			ev(j) = ev(j) + ( dcoef(k) - ccoef(k) ) * kv(j,k)
		enddo
		ev(j) = abs( tol / ev(j) ) ** ( 1.0 / real( order2 ) )
	enddo

	! Calculate the minimum error in the vector
	do j = 1,n
		if ( ev(j) < emin .OR. j == 1 ) then
			emin = ev(j)
		endif
	enddo

	! Use emin to calculate the new stepsize h
	s = upf * emin * h
	if ( s < hmax ) then
		h = s
	else
		h = hmax
	endif

	! Clear the k-vector kv, again
	do j = 1,n
		do k = 1,nstage
			kv(j,k) = 0.0
		enddo
	enddo

	! Perform the Runge-Kutta calculation again using new h
	do j = 1,nstage
		do k = 1,n
			xt(k) = x(k)
		enddo
		do k = 1,j
			do m = 1,n
				xt(m) = xt(m) + bcoef(k,j) * kv(m,k)
			enddo
		enddo
		tt = t + acoef(j) * h
		call grad (n,tt,xt,drv,data)
		do k = 1,n
			kv(k,j) = h * drv(k)
		enddo
	enddo
	do j = 1,n
		xt(j) = x(j) ! Copy then modify
		do k = 1,nstage
			xt(j) = x(j) + ccoef(k) * kv(j,k);
		enddo
	enddo
	t = t + h
end subroutine
