!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Numerical Recipes                                                    %
!                                                                      %
! William H. Press, et al, Numerical Recipes in Fortran 90 : The Art   %
! of Parallel Scientific Computing (Fortran Numerical Recipes, Vol 2), %
! Cambridge University Press, 1992, 571 pages.                         %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION splint(xa,ya,y2a,x,n)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY: locate
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
	REAL(SP), INTENT(IN) :: x
	INTEGER(I4B), OPTIONAL, INTENT(INOUT) :: n
	REAL(SP) :: splint
	INTEGER(I4B) :: khi,klo
	REAL(SP) :: a,b,h
    if (.not.present(n)) then
        n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    end if
	klo=max(min(locate(xa(1:n),x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
END FUNCTION splint
