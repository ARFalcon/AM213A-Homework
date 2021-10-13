program pi_appx

implicit none
integer, parameter :: dp = kind(1.d0)
real(kind=dp) :: appx_pi,true_pi,diff,n
integer :: i 
real(kind=dp), dimension(4) :: threshold 
threshold= (/1.e-4,1.e-8,1.e-12,1.e-16/)
n=0
appx_pi=0
diff=1.0
true_pi=acos(-1.d0)
do i=1,4
do while (diff>threshold(i))

    appx_pi=appx_pi+16**(-1*n)*(4/(8*n+1)-(2/(8*n+4))-(1/(8*n+5))-(1/(8*n+6)))
    diff=abs(appx_pi-true_pi)
    n=n+1

end do

write(*,*) 'After ',n,'iterations the approximate value of pi is',appx_pi,&
&'and the difference between the approximate value and true value is',diff
end do
end program pi_appx