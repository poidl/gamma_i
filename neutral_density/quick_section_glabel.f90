subroutine quick_section_glabel(s,t,p,longs,lats,nx,nz,g,glo,ghi)
	
implicit real*8(a-h,o-z)

data write_file/0/




integer :: nx, nz, ncast

real*8, dimension(nx) :: longs, lats

real*8, dimension(nz,nx) :: s, t, p, g, glo, ghi
	

g = -99.0d0
glo = -99.0d0
ghi = -99.0d0

if(write_file.eq.1) then
  open(11,file='section_dump.dat',status='unknown')
  write(11,*) nx,nz
  write(11,*) s,t,p,longs,lats
  close(11)
end if


do i = 1,nx
  if(s(1,i).GT.-99.0d0) then

    do k = 1,nz
      if(s(k,i).GT.-99.0d0) then
        ncast = k
      end if
    end do
    
!    open(11,file='data.dat',status='unknown')
!    write(11,*) nx,nz
!    write(11,*) longs(i),lats(i),ncast
!    write(11,'(3f12.4)') (s(k,i),t(k,i),p(k),k=1,ncast)
!    close(11)

    call gamma_n(s(1,i),t(1,i),p(1,i),ncast,longs(i),lats(i),g(1,i),glo(1,i),ghi(1,i))

  end if
end do


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
integer :: s_pr, t_pr, p_pr, longs_pr, lats_pr
integer :: g_pr, glo_pr, ghi_pr
!-----------------------------------------------------------------------
!

integer :: nlhs, nrhs, mxGetM, mxGetN
            
integer :: nx, nz, nxz, called

real*8, dimension(:),   allocatable :: longs, lats  

real*8, dimension(:,:), allocatable :: s, t, p, g, glo, ghi

save called
	
data called/0/


!    Check for proper number of arguments. 
if (nrhs .ne. 5) then
  call mexErrMsgTxt('5 inputs required.')
elseif (nlhs .ne. 3) then
  call mexErrMsgTxt('3 output required.')
endif
      
!    Get the nxz of the input matrices

nz = mxGetM(prhs(1))
nx = mxGetN(prhs(1))

nxz = nx*nz
      

!    allocat array sizes
      
if(called.eq.0) then
  allocate(s(1:nz,1:nx)); allocate(t(1:nz,1:nx))
  allocate(p(1:nz,1:nx)); allocate(g(1:nz,1:nx))
  allocate(glo(1:nz,1:nx)); allocate(ghi(1:nz,1:nx))
  allocate(longs(1:nx)); allocate(lats(1:nx))
  called = 1
endif

!    create matrices for the return arguments

plhs(1) = mxCreateFull(nz,nx,-99.0)
plhs(2) = mxCreateFull(nz,nx,-99.0)
plhs(3) = mxCreateFull(nz,nx,-99.0)

!    set up pointers

    s_pr  =  mxGetPr(prhs(1))
    t_pr  =  mxGetPr(prhs(2))
    p_pr  =  mxGetPr(prhs(3))
longs_pr  =  mxGetPr(prhs(4))
 lats_pr  =  mxGetPr(prhs(5))

    g_pr  =  mxGetPr(plhs(1))
  glo_pr  =  mxGetPr(plhs(2))
  ghi_pr  =  mxGetPr(plhs(3))


!    load the data into Fortran arrays

call mxCopyPtrToReal8(    s_pr,     s, nxz)
call mxCopyPtrToReal8(    t_pr,     t, nxz)
call mxCopyPtrToReal8(    p_pr,     p, nxz)
call mxCopyPtrToReal8(longs_pr, longs, nx)
call mxCopyPtrToReal8( lats_pr,  lats, nx)

	
!    call the computational subroutine

call quick_section_glabel(s,t,p,longs,lats,nx,nz,g,glo,ghi)
 	

!    load the output into a MATLAB array

call mxCopyReal8ToPtr(g,     g_pr, nxz)
call mxCopyReal8ToPtr(glo, glo_pr, nxz)
call mxCopyReal8ToPtr(ghi, ghi_pr, nxz)	


return
end     
