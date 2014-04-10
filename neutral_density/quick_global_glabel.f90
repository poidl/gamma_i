subroutine quick_global_glabel(s,t,p,longs,lats,nx,ny,nz,g,glo,ghi)
	
implicit real*8(a-h,o-z)

data write_file/0/

!	t is in-situ or conservative temperature depending on library



integer :: nx, ny, nz, ncast

real*8, dimension(ny, nx) :: longs, lats

 	
real*8, dimension(nz,ny,nx) :: s, t, p, g, glo, ghi
	

g = -99.0d0
glo = -99.0d0
ghi = -99.0d0

if(write_file.eq.1) then
  open(11,file='ocean_dump.dat',status='unknown')
  write(11,*) nx,ny,nz
  write(11,*) s,t,p,longs,lats
  close(11)
end if


do j = 1,ny
do i = 1,nx
  if(s(1,j,i).GT.-99.0d0.AND.-80.0d0.LE.lats(j,i).AND.lats(j,i).LE.90.0d0) then

    do k = 1,nz
      if(s(k,j,i).GT.-99.0d0) then
        ncast = k
      end if
    end do
    
!    open(11,file='data.dat',status='unknown')
!    write(11,*) nx,ny,nz
!    write(11,*) longs(i),lats(j),ncast
!    write(11,'(3f12.4)') (s(k,j,i),t(k,j,i),p(k),k=1,ncast)
!    close(11)

    call gamma_n(s(1,j,i),t(1,j,i),p(1,j,i),ncast,longs(j,i),lats(j,i),g(1,j,i),glo(1,j,i),ghi(1,j,i))

  end if
end do
end do	 



return
end




!     This is a MEX file for MATLAB.
!     Copyright (c) 1984-98 by The MathWorks, Inc.
!     $Revision: 1.7 $
      
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
	
implicit real*8(a-h,o-z)

!-----------------------------------------------------------------------
integer plhs(*), prhs(*)
integer :: mxCreateFull
integer :: s_pr, t_pr, p_pr, longs_pr, lats_pr
integer :: g_pr, glo_pr, ghi_pr
!-----------------------------------------------------------------------
!

integer :: nlhs, nrhs, mxGetM, mxGetN
            
integer :: nx, ny, nz, nxy, nxyz, called

real*8, dimension(:,:),   allocatable :: longs, lats  

real*8, dimension(:,:), allocatable :: s, t, p, g, glo, ghi

save called
	
data called/0/


!    Check for proper number of arguments. 
if (nrhs .ne. 5) then
  call mexErrMsgTxt('5 inputs required.')
elseif (nlhs .ne. 3) then
  call mexErrMsgTxt('3 output required.')
endif
      
!    Get the nxyz of the input matrices

nz = mxGetM(prhs(1))
nxy = mxGetN(prhs(1))
nxyz = nxy*nz
      
nx = mxGetM(prhs(4)); ny = mxGetM(prhs(5));


!    check to see dimensions balance

if (nxy .ne. nx*ny) then
  open(10,file='here.txt',status='unknown')
  write(10,*) nx,ny,nz
  write(10,*) nx*ny,nxy,nxyz
  close(10)
  call mexErrMsgTxt('Dimensions do not balance ... see here file')
endif


!    allocat array sizes
      
if(called.eq.0) then
  allocate(s(1:nz,1:nxy)); allocate(t(1:nz,1:nxy))
  allocate(p(1:nz,1:nxy)); allocate(g(1:nz,1:nxy))
  allocate(glo(1:nz,1:nxy)); allocate(ghi(1:nz,1:nxy))
  allocate(longs(1:ny,1:nx)); allocate(lats(1:ny,1:nx))
  called = 1
endif

!    create matrices for the return arguments

plhs(1) = mxCreateFull(nz,nxy,0)
plhs(2) = mxCreateFull(nz,nxy,0)
plhs(3) = mxCreateFull(nz,nxy,0)

!    set up pointers

s_pr =  mxGetPr(prhs(1))
t_pr =  mxGetPr(prhs(2))
p_pr =  mxGetPr(prhs(3))
longs_pr =  mxGetPr(prhs(4))
lats_pr =  mxGetPr(prhs(5))

g_pr = mxGetPr(plhs(1))
glo_pr = mxGetPr(plhs(2))
ghi_pr = mxGetPr(plhs(3))


!    load the data into Fortran arrays

call mxCopyPtrToReal8(s_pr,s,nxyz)
call mxCopyPtrToReal8(t_pr,t,nxyz)
call mxCopyPtrToReal8(p_pr,p,nxyz)
call mxCopyPtrToReal8(longs_pr,longs,nxy)
call mxCopyPtrToReal8(lats_pr,lats,nxy)

	
!    call the computational subroutine

call quick_global_glabel(s,t,p,longs,lats,nx,ny,nz,g,glo,ghi)
 	

!    load the output into a MATLAB array

call mxCopyReal8ToPtr(g,g_pr,nxyz)
call mxCopyReal8ToPtr(glo,glo_pr,nxyz)
call mxCopyReal8ToPtr(ghi,ghi_pr,nxyz)	


return
end     
