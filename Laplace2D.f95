program Laplace2D

implicit none

	integer :: Nv, Nvf, Ne, a, b, c, Nvv, i, j
	real :: D, total, add, saa, sbb, scc, sab, sbc, sca, D1, D2, sab1, sab2, sbc1,sbc2,sca1,sca2
	real, dimension(:), allocatable :: U, ConstantVector, UVector
	real, dimension(:,:), allocatable :: location, S
	integer, dimension(:,:), allocatable :: element
	real, dimension(:,:), allocatable :: CoefficientMatrix, GaussianMatrix

	interface 
      		subroutine SolveLinearEquation (M)
         		real, dimension(:,:), intent (inout) :: M       
      		end subroutine SolveLinearEquation      
  	end interface 

	!open a file
	open(unit = 1, file = "input.txt", form = "formatted", status = "old", action = "read")

	!read the number of boundary vertices
	read (unit  = 1, fmt  = *) Nvf

	!read the total number of vertices
	read (unit  = 1, fmt  = *) Nv

	allocate ( U (Nv) )
	allocate (location (Nv,2))

	!read the coordinates of the vertices
	do i = 1, Nv
		read (unit  = 1, fmt  = *) location(i,1), location(i,2)	
	end do

	!read the values of U at the boundary vertices
	do i = 1, Nvf
		read (unit  = 1, fmt  = *) U(i)
	end do

	!read the number of the triangular elements in the mesh
	read (unit  = 1, fmt  = *) Ne

	!create a Nex3 matrix called element which stores the information about the triangles
	allocate ( element (Ne,3) )

	!read the triangular elements
	do i = 1, Ne
		read (unit  = 1, fmt  = *) element(i,1), element(i,2), element(i,3)
	end do

	allocate ( S(Nv,Nv) )	
	S(:,:) = 0
	
	!traverse throught the elements and add the corresponding entries into S matrix
	do i = 1, Ne
		a = element(i,1)
		b = element(i,2)
		c = element(i,3)
		D1 = location(b,1)*location(c,2) - location(c,1)*location(b,2) 
		D2 = location(c,1)*location(a,2) - location(a,1)*location(c,2) + location(a,1)*location(b,2) - location(b,1)*location(a,2)
		D = D1 + D2
		saa = ( (location(b,2) - location(c,2))**2 + (location(c,1) - location(b,1))**2 ) / D
		sbb = ( (location(c,2) - location(a,2))**2 + (location(a,1) - location(c,1))**2 ) / D
		scc = ( (location(a,2) - location(b,2))**2 + (location(b,1) - location(a,1))**2 ) / D
		sab1 = (location(b,2) - location(c,2))*(location(c,2) - location(a,2))
		sab2 = (location(c,1) - location(b,1))*(location(a,1) - location(c,1))
		sab = (sab1 + sab2)/D
		sbc1 = (location(c,2) - location(a,2))*(location(a,2) - location(b,2))
		sbc2 = (location(a,1) - location(c,1))*(location(b,1) - location(a,1))
		sbc = (sbc1 + sbc2)/D
		sca1 = (location(a,2) - location(c,2))*(location(b,2) - location(c,2))
		sca2 = (location(c,1) - location(a,1))*(location(c,1) - location(b,1))
		sca = (sca1 + sca2)/D
		S(a,a) = S(a,a) + saa;
		S(b,b) = S(b,b) + sbb;
		S(c,c) = S(c,c) + scc;
		S(a,b) = S(a,b) + sab;
		S(b,c) = S(b,c) + sbc;
		S(c,a) = S(c,a) + sca;	
	end do
	
	Nvv = Nv - Nvf
	allocate (GaussianMatrix(Nvv,Nvv+1))
	allocate (CoefficientMatrix(Nvv,Nvv))
	allocate (ConstantVector(Nvv))
	allocate (UVector(Nvv))

	!calculate CoefficientMatrix and ConstantVector
	do i = Nvf+1, Nv
		total = 0
		do j  = 1, Nvf
			add = S(j,i) + S(i,j)
			total = total + add*U(j)
		end do
		ConstantVector(i-Nvf) = -1 * total
		do j = Nvf+1, Nv
			CoefficientMatrix(i-Nvf,j-Nvf) = S(i,j) + S(j,i)
		end do
	end do 
 
	!Solve linear eqautions using Gauss Jordan Method
	GaussianMatrix(:,1:Nvv) = CoefficientMatrix
	GaussianMatrix(:,Nvv+1) = ConstantVector

	call SolveLinearEquation(GaussianMatrix)

	Uvector = GaussianMatrix(:,Nvv+1)
	do i = 1, Nvv
		print *, Uvector(i)
	end do

end program Laplace2D


subroutine SolveLinearEquation(M)

implicit none

	real, dimension(:,:), intent(inout) :: M
	integer :: n, i, j, s(2)
	real :: c, b
	s = shape(M)
	n = s(1)
	
	do i = 1, n

		b = M(i,i)
		M(i,:) = M(i, :)/b

		do j = 1, n
			if (j /= i)	then
				if(M(j,i) == 0)		then
					cycle
				else
					b = M(j,i)
					M(j,:) = M(j,: ) - M(i,: ) *b
				endif
			endif
		end do
	end do

end subroutine SolveLinearEquation 
