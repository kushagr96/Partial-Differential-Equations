program Laplace1D

implicit none

	integer :: Nv, Ne, i, ok
	integer, dimension(:), allocatable :: pivot
	real, dimension(:), allocatable :: NodesPosition, ConstantVector, UVector
	real, dimension(:,:), allocatable :: CoefficientMatrix, GaussianMatrix
	real :: width, Ua, Ub, a, b
	

	interface 
      		subroutine SolveLinearEquation (M)
         		real, dimension(:,:), intent (inout) :: M       
      		end subroutine SolveLinearEquation      
  	end interface 

	print *, "Enter the domain: "
	read *, a, b

	print *, "Enter the number of nodes excluding the end points: "
	read *, Nv
	
	
	print *, "Enter the boundary conditions: "
	read *, Ua, Ub

	allocate ( NodesPosition(Nv+2) )
	Ne = Nv + 1
	width  = (b-a)/Ne
	
	do i = 1, Ne
		NodesPosition(i) = a + width*(i-1)
	end do
	NodesPosition(Nv+2) = b
	
	print *, "The position of the nodes is:"
	do i = 1, Ne+1
		print *, NodesPosition(i)
	end do

	allocate( CoefficientMatrix(Nv, Nv) )
	allocate (ConstantVector(Nv))
	allocate ( GaussianMatrix(Nv,Nv+1))
	allocate (pivot(Nv))
	
	ConstantVector(:) = 0
	ConstantVector(1) = ConstantVector(1) + Ua
	ConstantVector(Nv) = ConstantVector(Nv) + Ub
	CoefficientMatrix(:,:) = 0
	
	do i = 1, Nv
		CoefficientMatrix(i, i) = 2
	end do

	do i = 1, Nv
		if ( i-1 > 0 )	then
			CoefficientMatrix(i, i-1) = -1
		endif
		if ( i+1 <= Nv )  then
			CoefficientMatrix(i, i+1) = -1	
		endif
	end do

	GaussianMatrix(:,1:Nv) = CoefficientMatrix
	GaussianMatrix(:,Nv+1) = ConstantVector
	!call SGESV(Nv, 1, CoefficientMatrix, Nv, pivot, ConstantVector, Nv, ok)

	call SolveLinearEquation(GaussianMatrix)

	Uvector = GaussianMatrix(:,Nv+1)

	print *, "The values of the function at the nodes is: "	
	do i = 1, Nv
		print *, i, ":", UVector(i) 
	end do	
	
end program Laplace1D



subroutine SolveLinearEquation(M)

implicit none

	real, dimension(:,:), intent(inout) :: M
	integer :: n, i, j, s(2)
	real :: c, b
	real, dimension(:), allocatable :: v, v2
	s = shape(M)
	n = s(1)

	allocate (v(n) )
	allocate (v2(n) )

	do i = 1, n

		b = M(i,i)
		v = M(i, :)/b
		M(i,:) = v

		do j = 1, n
			if (j /= i)	then
				if(M(j,i) == 0)		then
					cycle
				else
					b = M(j,i)
					v = M(i,:) * b
					v2 = M(j,:) - v
					M(j,:) = v2 
				endif
			endif
		end do
	end do

end subroutine SolveLinearEquation 
