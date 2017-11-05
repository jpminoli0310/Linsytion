!---------------------------Header-------------------------------!
! MODULE Lynsition
! Libreria que contiene subrutinas para solución de sistemas de
! ecuaciones lineales.
! 
! AUTOR: Juan Pablo Velásquez Minoli  
!-------------------------End Header-----------------------------!

MODULE Lynsition
IMPLICIT NONE
!        ========
CONTAINS
!        ========
	SUBROUTINE conjugateGradient(A,b,cond,x,N,tol,nx)
		!----------------Input Variables----------------!		
		INTEGER                   , INTENT(IN)   :: N,nx
		REAL, DIMENSION(1:nx,1:nx), INTENT(IN)   :: A, cond
		REAL, DIMENSION(1:nx)     , INTENT(IN)   :: b
		REAL                      , INTENT(IN)   :: tol
		!----------------Work Variables----------------!
		INTEGER                   		 :: i,j
		REAL, DIMENSION(:)  , ALLOCATABLE        :: r,w,v,u
		REAL                                     :: alpha,t,beta,s
		REAL, DIMENSION(:)              	 :: x
		!----------------Init--------------------------!
		ALLOCATE(r(1:nx))
		ALLOCATE(w(1:nx))
		ALLOCATE(v(1:nx))
		ALLOCATE(u(1:nx))
		
		r=b-MATMUL(A,x)
		w=MATMUL(cond,r)
		v=MATMUL(cond,w)
		alpha=DOT_PRODUCT(w,w)

		PRINT *,"      Iteración      Solución"
		DO i=1,N
			IF (SQRT(DOT_PRODUCT(v,v))<=tol) THEN
				PRINT *, "Se ha encontrado una solución"
				PRINT *, "La solución es: ",x
				EXIT
			END IF
			u=MATMUL(A,v)
			t=alpha/DOT_PRODUCT(v,u)
			x=x+(t*v)
			PRINT *, i,"      ",x
			r=r-(t*u)
			w=MATMUL(cond,r)
			beta=DOT_PRODUCT(w,w)
			IF (SQRT(beta*beta)<=tol) THEN
				IF (SQRT(DOT_PRODUCT(r,r))<=tol) THEN
					PRINT *, "Se ha encontrado una solución"
					PRINT *, "La solución es: ",x
					EXIT
				END IF
			END IF
			s=beta/alpha
			v=(MATMUL(cond,w))+(s*v)
			alpha=beta
		END DO




		IF(ALLOCATED( r ) )    DEALLOCATE( r )
		IF(ALLOCATED( w ) )    DEALLOCATE( w )         
		IF(ALLOCATED( v ) )    DEALLOCATE( v )
		IF(ALLOCATED( u ) )    DEALLOCATE( u )
		
	END SUBROUTINE conjugateGradient
END MODULE Lynsition
