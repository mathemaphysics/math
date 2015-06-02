!> \brief Performs matrix multiplication

!> This function is a custom version of a LAPACK routine
!> which multiplies two matrices. To summarize,
!> \c C := \c alpha \c * \c A*B + \c beta \c * \c C
!> \param tA Transposition of \c A
!> \param tB Transposition of \c B
!> \param n Number rows in A
!> \param k Number of columns in \c A and rows in \c B
!> \param m Number of columns in \c B
!> \param alpha Multiplier on product \c AB
!> \param A Pointer to the column major matrix \c A
!> \param lda The leading dimension of \c A
!> \param B Pointer to the column major matrix \c B
!> \param ldb The leading dimension of \c B
!> \param beta Multiplier on additive \c C
!> \param C Pointer to the column major matrix \c C
!> \param ldc The leading dimension of \c C
subroutine dgemm(tA,tB,n,k,m,alpha,A,lda,B,ldb,beta,C,ldc)
        implicit none
        integer :: n,k,m,lda,ldb,ldc,i,j,p
	real*8, dimension(lda,n) :: A
	real*8, dimension(ldb,k) :: B
	real*8, dimension(ldc,n) :: C
	real*8 alpha,beta
        character :: tA,tB
        if ( tA == 't' .or. tA == 'T' ) then
                if( tB == 't' .or. tB == 'T' ) then
                        do i = 1,n
				do j = 1,m
					C(i,j) = beta * C(i,j)
					do p = 1,k
						C(i,j) = C(i,j) + alpha * A(p,i) * B(j,p)
					enddo
				enddo
			enddo
		else
			do i = 1,n
				do j = 1,m
					C(i,j) = beta * C(i,j)
					do p = 1,k
						C(i,j) = C(i,j) + alpha * A(p,i) * B(p,j)
					enddo
				enddo
			enddo
		endif
	else
		if( tB == 't' .or. tB == 'T' ) then
			do i = 1,n
				do j = 1,m
					C(i,j) = beta * C(i,j)
					do p = 1,k
						C(i,j) = C(i,j) + alpha * A(i,p) * B(j,p)
					enddo
				enddo
			enddo
		else
			do i = 1,n
				do j = 1,m
					C(i,j) = beta * C(i,j)
					do p = 1,k
						C(i,j) = C(i,j) + alpha * A(i,p) * B(p,j)
					enddo
				enddo
			enddo
		endif
	endif
end subroutine

subroutine dgemv(tA,n,m,alpha,A,lda,x,ldx,beta,y,ldy)
	implicit none
	integer :: i,k,n,m,lda,ldx,ldy
	real*8, dimension(lda,m) :: A
	real*8, dimension(ldx,m) :: x
	real*8, dimension(ldy,n) :: y
	real*8 :: alpha,beta
	character :: tA
	if( tA == 't' .or. tA == 'T' ) then
		do i = 1,n
			do k = 1,m
				y(1,i) = y(1,i) + A(i,k) * x(1,k)
			enddo
		enddo
	else
		do i = 1,n
			do k = 1,m
				y(1,i) = y(1,i) + A(k,i) * x(1,k)
			enddo
		enddo
	endif
end subroutine

subroutine daxpy(n,x,ldx,y,ldy,ac)
	implicit none
	integer :: n,ldx,ldy,i
	real*8, dimension(ldx,n) :: x ! Ensure that each entry of x is ldx quadwords apart
	real*8, dimension(ldy,n) :: y ! Ensure that each entry of y is ldy quadwords apart
	real*8 :: ac
	ac = 0.0
	do i = 1,n
		ac = ac + x(1,i) * y(1,i)
	enddo
end subroutine

subroutine dlacpy(uplo,n,m,A,lda,B,ldb)
	integer :: i,j,n,m,lda,ldb
	real*8, dimension(lda,m) :: A
	real*8, dimension(ldb,m) :: B
	character :: uplo
	if( uplo == 'u' .or. uplo == 'U' ) then
		do i = 1,m
			do j = 1,min(i,n)
				B(j,i) = A(j,i)
			enddo
		enddo
	else
		if ( uplo == 'l' .or. uplo == 'L' ) then
			do i = 1,n
				do j = 1,min(i,m)
					B(i,j) = A(i,j)
				enddo
			enddo
		else
			do i = 1,n
				do j = 1,m
					B(i,j) = A(i,j)
				enddo
			enddo
		endif
	endif
end subroutine

!> \brief Determines the permutation of rows/columns required for all nonzero
!> diagonals

!> This function determines the permutation matrix, \c P , such
!> that \c PA has no zeroes on the diagonal. The variable \c ipiv
!> is a one-dimensional vector mapping row/column numbers to their
!> new locations. This is primarily used internally to deal with
!> the LU factorization which requires that all diagonal elements
!> be nonzero.
!> \param n The dimension of the matrix \c A
!> \param A The square matrix of interest
!> \param lda The leading dimension of \c A
!> \param ipiv The pointer to the length \c n vector to store the permutatiion
!> \param res The result of the operation; returns 0 if all is well
subroutine dgepiv(n,A,lda,ipiv,res)
	implicit none
	integer :: n,lda,res,i,j,found,tmp
	integer, dimension(n) :: ipiv
	real*8, dimension(lda,n) :: A
	real*8, parameter :: zlim = 1.0e-6
	res = 0
	do i = 1,n ! Initialize all rows to their current locations
		ipiv(i) = i
	enddo
	do i = 1,n
		if( abs(A(ipiv(i),i)) < zlim ) then
			found = 0
			do j = 1,n
				if( ipiv(j) == ipiv(i) ) then
					continue
				else
					if( abs(A(ipiv(i),j)) > zlim .and. abs(A(ipiv(j),i)) > zlim ) then
						tmp = ipiv(i)
						ipiv(i) = ipiv(j)
						ipiv(j) = tmp
						found = 1
						exit
					endif
				endif
			enddo
			if( found == 0 ) then ! If there is even one zero which cannot be eliminated then matrix is rank deficient
				res = -1
				exit
			else
				res = res + 1 ! Count the number of zeros originally on the diagonal
			endif
		endif
	enddo
end subroutine

subroutine dgelu(perm,n,m,A,lda,ipiv,res)
	implicit none
	integer :: n,m,lda,res,i,j,k
	integer, dimension(n) :: ipiv
	real*8, dimension(lda,n) :: A
	real*8 :: tau
	character :: perm
	if( perm == 'n' .or. perm == 'N' ) then	! perm = 'n' or 'N' indicates DO NOT calculate the row permutation
		do i = 1,min(n,m)		! matrix, P, and also do not use the input ipiv; ipiv can just be NULL
			do j = i+1,n
				tau = A(j,i) / A(i,i)
				do k = i+1,m
					A(j,k) = A(j,k) - tau * A(i,k)
				enddo
				A(j,i) = tau
			enddo
		enddo
	else
		if( perm == 'f' .or. perm == 'F' ) then ! Call dgepiv to calculate the pivots and put it in ipiv as output
			call dgepiv(n,A,lda,ipiv,res)	! This function deals with square matrices, but it works for rectangular as well (I hope)
			if( res < 0 ) then		! This would indicate that the pivot matrix calculation failed
				return
			endif
			do i = 1,min(n,m)
				do j = i+1,n
					tau = A(ipiv(j),i) / A(ipiv(i),i)
					do k = i+1,m
						A(ipiv(j),k) = A(ipiv(j),k) - tau * A(ipiv(i),k)
					enddo
					A(ipiv(j),i) = tau
				enddo
			enddo
		else
			if( perm == 'p' .or. perm == 'P' ) then	! Use the ipiv pivots given as input
				do i = 1,min(n,m)
					do j = i+1,n
						tau = A(ipiv(j),i) / A(ipiv(i),i)
						do k = i+1,m
							A(ipiv(j),k) = A(ipiv(j),k) - tau * A(ipiv(i),k)
						enddo
						A(ipiv(j),i) = tau
					enddo
				enddo
			endif
		endif
	endif
end subroutine

subroutine dgefb(perm,n,nrhs,A,lda,ipiv,B,ldb,res)
	implicit none
	integer :: n,nrhs,lda,ldb,res,i,j,k
	integer, dimension(n) :: ipiv
	real*8, dimension(lda,n) :: A
	real*8, dimension(ldb,nrhs) :: B
	character :: perm
	if( perm == 'n' .or. perm == 'N' ) then
		do k = 1,nrhs
			do i = 2,n ! Forward substitution
				do j = 1,i-1
					B(i,k) = B(i,k) - A(i,j) * B(j,k)
				enddo
			enddo
			do i = n,1,-1 ! Backward substitution
				do j = n,i+1,-1
					B(i,k) = B(i,k) - A(i,j) * B(j,k)
				enddo
				B(i,k) = B(i,k) / A(i,i)
			enddo
		enddo
		res = 0
	else
		if( perm == 'f' .or. perm == 'F' ) then
			call dgepiv(n,A,lda,ipiv,res)
			if( res >= 0 ) then
				do k = 1,nrhs
					do i = 2,n ! Forward substitution
						do j = 1,i-1
							B(ipiv(i),k) = B(ipiv(i),k) - A(ipiv(i),j) * B(ipiv(j),k)
						enddo
					enddo
					do i = n,1,-1 ! Backward substitution
						do j = n,i+1,-1
							B(ipiv(i),k) = B(ipiv(i),k) - A(ipiv(i),j) * B(ipiv(j),k)
						enddo
						B(ipiv(i),k) = B(ipiv(i),k) / A(ipiv(i),i)
					enddo
				enddo
			endif
			res = 0
		else
			if( perm == 'p' .or. perm == 'P' ) then
				do k = 1,nrhs
					do i = 2,n ! Forward substitution
						do j = 1,i-1
							B(ipiv(i),k) = B(ipiv(i),k) - A(ipiv(i),j) * B(ipiv(j),k)
						enddo
					enddo
					do i = n,1,-1 ! Backward substitution
						do j = n,i+1,-1
							B(ipiv(i),k) = B(ipiv(i),k) - A(ipiv(i),j) * B(ipiv(j),k)
						enddo
						B(ipiv(i),k) = B(ipiv(i),k) / A(ipiv(i),i)
					enddo
				enddo
				res = 0
			endif
		endif
	endif
end subroutine

subroutine dgedet(perm,n,A,lda,ipiv,det,res)
	implicit none
	integer :: n,lda,res,i,j,tmp
	integer, dimension(n) :: ipiv
	real*8, dimension(lda,n) :: A
	real*8 :: det
	character :: perm
	if( n > 0 ) then
		if( n == 1 ) then
			det = A(1,1)
			res = 0
		else
			res = 0
			if( perm == 'f' .or. perm == 'F' ) then
				call dgepiv(n,A,lda,ipiv,res)
			endif
			if( res >= 0 ) then
				call dgelu(perm,n,n,A,lda,ipiv,res)
				if( res >= 0 ) then
					det = 1.0
					if( perm .ne. 'n' .and. perm .ne. 'N' ) then
						do i = 1,n
							det = det * A(ipiv(i),i)
						enddo
						tmp = 0
						do i = 1,n
							do j = i+1,n
								if( ipiv(i) > ipiv(j) ) then
									tmp = tmp + 1
								endif
							enddo
						enddo
						tmp = mod(tmp,2)
						if( tmp == 1 ) then
							det = -1.0 * det
						endif
					else
						do i = 1,n
							det = det * A(i,i)
						enddo
					endif
				endif
			endif
		endif
	else
		res = -1
	endif
end subroutine

!> \brief Solves the linear system \c AX = \c B, overwriting \c B

!> This function is analogous to the LAPACK routine
!> by the same name. It solves the system \c AX \c = \c B
!> using LU factorization. The variable \c ipiv is space
!> required to calculate the necessary permutation of rows
!> and columns to assure that no diagonal entries are zero.
!> The solution is stored in \c B, overwriting the original
!> input.
!> \param n The number of rows (and columns) in \c A
!> \param nrhs The number of columns in \c B
!> \param A The matrix to invert
!> \param lda The leading dimension of \c A
!> \param ipiv The permutation vector
!> \param B The righthand side of the system
!> \param ldb The leading dimension of \c B
!> \param res The result code; returns 0 for all okay
subroutine dgesv(n,nrhs,A,lda,ipiv,B,ldb,res)
	implicit none
	integer :: n,nrhs,lda,ldb,ldx,res,k,i
	integer, dimension(n) :: ipiv
	real*8, dimension(n) :: temp
	real*8, dimension(lda,n) :: A
	real*8, dimension(ldb,nrhs) :: B
	call dgelu('f',n,n,A,lda,ipiv,res)
	if( res >= 0 ) then
		call dgefb('p',n,nrhs,A,lda,ipiv,B,ldb,res)
		if( res >= 0 ) then
			do k = 1,nrhs		! Putting the values into the right places in each row
				do i = 1,n	! Row i is physically located at row ipiv(i) in memory
					temp(i) = B(i,k)
				enddo
				do i = 1,n
					B(i,k) = temp(ipiv(i))
				enddo
			enddo
		endif
	endif
end subroutine

subroutine matrix_print(n,m,A)
	implicit none
	real*8, dimension(n,m) :: A
	integer :: n,m
	real*8 x
	integer i,j
	do i = 1,n
		do j = 1,m
			x = A(i,j)
			write (*,"(2f11.7)",advance='no') x
		enddo
		print *
	enddo
end subroutine
