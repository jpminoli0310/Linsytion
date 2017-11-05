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
		REAL, DIMENSION(1:nx)     , INTENT(IN)   :: b, x
		REAL                      , INTENT(IN)   :: tol
		!----------------Work Variables----------------!
		INTEGER                   		 :: i,j

	END SUBROUTINE conjugateGradient
END MODULE Lynsition
