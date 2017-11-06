!----------------------Example1---------------------------------!
!Programa ejemplo para el uso de la subrutina conjugateGradient
!para la solucion de un sistema de ecuaciones lineales de la
!la forma Ax=b
!
!USES   : Linsytion
!---------------------------------------------------------------!
PROGRAM Example1
USE Linsytion   
IMPLICIT NONE
REAL, DIMENSION(:,:), ALLOCATABLE  :: A, cond
REAL, DIMENSION(:)  , ALLOCATABLE  :: b, x
REAL                               :: tol
INTEGER                            :: N,nx
INTEGER                            :: startTime
INTEGER                            :: endTime      
INTEGER                            :: countRate
REAL*8                             :: elapsedTime

CALL SYSTEM_CLOCK(startTime, countRate) ! Inicia Reloj
N=3
tol=0.05
nx=3

ALLOCATE(cond(1:nx,1:nx))
ALLOCATE(A(1:nx,1:nx))
ALLOCATE(b(1:nx))
ALLOCATE(x(1:nx))

A     = RESHAPE((/ 4, 3, 0, 3, 4, -1, 0, -1, 4 /),SHAPE=(/nx,nx/))
cond  = RESHAPE((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /),SHAPE=(/nx,nx/))
b     = (/ 24, 30, -24 /)
x     = 0

call conjugateGradientPara(A,b,cond,x,N,tol,nx)
PRINT *, "La solución al sistema es: ",x

IF(ALLOCATED( A ) )     DEALLOCATE( A )
IF(ALLOCATED( cond ) )    DEALLOCATE( cond )         
IF(ALLOCATED( b ) )     DEALLOCATE( b )         
IF(ALLOCATED( x  ) )    DEALLOCATE( x  )

CALL SYSTEM_CLOCK(endTime)              ! Finaliza Reloj
ElapsedTime= REAL(endTime - startTime) / REAL(countRate)

PRINT *, '******************************************************'
PRINT *, 'TIME', elapsedTime
PRINT *, 'Done.'
END PROGRAM
