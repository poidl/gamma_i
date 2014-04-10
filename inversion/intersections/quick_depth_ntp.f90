subroutine quick_depth_ntp(s0,t0,p0,s,t,p,n,eos,sntp,tntp,pntp)
	
implicit real*8(a-h,o-z)

data write_file/0/

!	t is in-situ, potential or conservative temperature depending on eos


integer :: n

real*8 :: s0, t0, p0, eos, sntp, tntp, pntp

real*8, dimension(n) :: s, t, p

	

call depth_ns(s,t,p,n,s0,t0,p0,sntp,tntp,pntp)


return
end




!     This is a MEX file for MATLAB.
!     Copyright (c) 1984-98 by The MathWorks, Inc.
!     $Revision: 1.7 $
      
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	
implicit real*8(a-h,o-z)

!-----------------------------------------------------------------------
integer    plhs(*), prhs(*)
integer :: mxCreateFull
integer :: s0_pr, t0_pr, p0_pr
integer :: s_pr, t_pr, p_pr, eos_pr
integer :: sntp_pr, tntp_pr, pntp_pr
!-----------------------------------------------------------------------
!

integer :: nlhs, nrhs, mxGetM, mxGetN
            
integer :: n

real*8 :: s0, t0, p0, eos, sntp, tntp, pntp

real*8, dimension(:),   allocatable :: s, t, p


!    Check for proper number of arguments. 
if (nrhs .ne. 7) then
  call mexErrMsgTxt('7 inputs required.')
elseif (nlhs .ne. 3) then
  call mexErrMsgTxt('3 output required.')
endif
      
!    Get the n of the input vectors

n = mxGetM(prhs(4)); n1 = 1
!nx = mxGetN(prhs(1))
      

!    allocat array sizes
      
allocate(s(1:n)); allocate(t(1:n)); allocate(p(1:n))


!    create matrices for the return arguments

plhs(1) = mxCreateFull(n1,n1,-99.0)
plhs(2) = mxCreateFull(n1,n1,-99.0)
plhs(3) = mxCreateFull(n1,n1,-99.0)

!    set up pointers

 s0_pr  =  mxGetPr(prhs(1))
 t0_pr  =  mxGetPr(prhs(2))
 p0_pr  =  mxGetPr(prhs(3))
  s_pr  =  mxGetPr(prhs(4))
  t_pr  =  mxGetPr(prhs(5))
  p_pr  =  mxGetPr(prhs(6))
eos_pr  =  mxGetPr(prhs(7))

sntp_pr  =  mxGetPr(plhs(1))
tntp_pr  =  mxGetPr(plhs(2))
pntp_pr  =  mxGetPr(plhs(3))


!    load the data into Fortran arrays

call mxCopyPtrToReal8( s0_pr,  s0, n1)
call mxCopyPtrToReal8( t0_pr,  t0, n1)
call mxCopyPtrToReal8( p0_pr,  p0, n1)
call mxCopyPtrToReal8(  s_pr,   s, n)
call mxCopyPtrToReal8(  t_pr,   t, n)
call mxCopyPtrToReal8(  p_pr,   p, n)
call mxCopyPtrToReal8(eos_pr, eos, n1)

	
!    call the computational subroutine

call quick_depth_ntp(s0,t0,p0,s,t,p,n,eos,sntp,tntp,pntp) 	

!    load the output into a MATLAB array

call mxCopyReal8ToPtr(sntp, sntp_pr, n1)
call mxCopyReal8ToPtr(tntp, tntp_pr, n1)
call mxCopyReal8ToPtr(pntp, pntp_pr, n1)	


deallocate(s,t,p,stat=ierr)


return
end     
