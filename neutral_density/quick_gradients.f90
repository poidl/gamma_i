subroutine quick_gradients(s,t,p,g,ocean,en,longs,lats,opts,nx,ny,nz,nsx,nsy)


!	NOTE: t is in-situ temperature, depending on libraries

	
implicit real*8(a-h,o-z)

integer :: nx, ny, nz, write_file

integer, dimension(:,:), allocatable :: iocean, n

real*8, dimension(2) :: opts

!real*8, dimension(nx) :: longs
!real*8, dimension(ny) :: lats

!real*8, dimension(:,:), allocatable :: longs2, lats2 

real*8, dimension(ny,nx) :: longs, lats, ocean, en
 	
real*8, dimension(nz,ny,nx) :: s, t, p, g, nsx, nsy

common /ifd/ icalled_from_driver
	
save called
	
data called/0/, write_file/0/

  
icalled_from_driver = 0

!    allocate array sizes
      
if(called.eq.0) then
!  allocate(longs2(1:ny,1:nx)); allocate(lats2(1:ny,1:nx))
  allocate(iocean(1:ny,1:nx)); allocate(n(1:ny,1:nx))
  called = 1
endif	


iocean = int(ocean); n = int(en); iopt = int(opts(1))

if(write_file.eq.1) then
    open(11,file='ocean_dump.dat',status='unknown')
    write(11,*) nx,ny,nz
    write(11,*) s,t,p,g,iocean,longs,lats
    close(11)
end if

if(iopt.eq.0) then
  open(11,file='data/ocean_dump.dat',status='unknown')
  write(11,'(i4)') nx,ny,nz
  write(11,'(e20.12)') longs,lats,p,s,t,g,ocean
  close(11)
elseif(iopt.eq.1) then
!  do j = 1,ny
!  do i = 1,nx
!	longs2(j,i) = longs(i); lats2(j,i) = lats(j)
!  end do
!  end do
  call sigl_gradients(s,t,p,ocean,longs,lats,nx,ny,nz,nsx,nsy);
elseif(iopt.eq.2) then
!  do j = 1,ny
!  do i = 1,nx
!	longs2(j,i) = longs(i); lats2(j,i) = lats(j)
!  end do
!  end do
  call ns_gradients(s,t,p,g,ocean,longs,lats,nx,ny,nz,nsx,nsy)
elseif(iopt.eq.3) then
  do j = 1,ny
  do i = 1,nx
	if(p(nz,j,i).gt.0.) then
		i0 = i; j0 = j
	end if
  end do
  end do
!  call gpoly_gradients(s,t,p(1,j0,i0),longs,lats,iocean,nx,ny,nz,nsx,nsy)
elseif(iopt.eq.4) then
!  do j = 1,ny
!  do i = 1,nx
!	if(p(nz,j,i).gt.0.) then
!		i0 = i; j0 = j
!	end if
!  end do
!  end do
  pr = opts(2)
  call sigp_gradients(s,t,p,longs,lats,pr,nx,ny,nz,nsx,nsy)
elseif(iopt.eq.5) then
!  call ortho_gradients(s,t,p,longs,lats,nx,ny,nz,nsx,nsy)
end if


return
end


!     This is a MEX file for MATLAB
      
subroutine mexFunction(nlhs, plhs, nrhs, prhs)


integer plhs(nlhs), prhs(nrhs)
integer :: mxCreateFull
integer :: s_pr, t_pr, p_pr, g_pr, ocean_pr, n_pr, longs_pr, lats_pr, opts_pr
integer :: nsx_pr, nsy_pr


integer nlhs, nrhs
integer mxGetM, mxGetN
            
integer :: nx, ny, nz, nxy, nxyz

real*8, dimension(2) :: opts 
 	 
real*8, dimension(:,:), allocatable :: longs, lats  

real*8, dimension(:,:), allocatable :: s, t, p, g, ocean, n, nsx, nsy

save called
	
data called/0/


!    check for proper number of arguments

if (nrhs .ne. 9) then
  call mexErrMsgTxt('9 inputs required.')
elseif (nlhs .ne. 2) then
  call mexErrMsgTxt('2 output required.')
endif

      
!    get the nxyz of the input matrices

nz = mxGetM(prhs(1))
nxy = mxGetN(prhs(1))
nxyz = nxy*nz
      
nx = mxGetN(prhs(5)); ny = mxGetM(prhs(5))


!    check to see dimensions balance

if(nxy.ne.nx*ny) then
  open(10,file='here.txt',status='unknown')
  write(10,*) nx,ny,nz
  write(10,*) nx*ny,nxy,nxyz
  close(10)
  call mexErrMsgTxt('Dimensions do not balance ... see here file')
endif


!    allocate array sizes
      
if(called.eq.0) then
  allocate(longs(1:ny,1:nx)); allocate(lats(1:ny,1:nx));
  allocate(s(1:nz,1:nxy)); allocate(t(1:nz,1:nxy))
  allocate(p(1:nz,1:nxy)); allocate(g(1:nz,1:nxy))
  allocate(ocean(1:ny,1:nx)); allocate(n(1:ny,1:nx))
  allocate(nsx(1:nz,1:nxy)); allocate(nsy(1:nz,1:nxy))
  called = 1
endif


!    create matrices for the return arguments

plhs(1) = mxCreateFull(nz,nxy,0)
plhs(2) = mxCreateFull(nz,nxy,0)


!    set up pointers

s_pr =  mxGetPr(prhs(1))
t_pr =  mxGetPr(prhs(2))
p_pr =  mxGetPr(prhs(3))
g_pr = mxGetPr(prhs(4))
ocean_pr = mxGetPr(prhs(5))
n_pr = mxGetPr(prhs(6))
longs_pr =  mxGetPr(prhs(7))
lats_pr =  mxGetPr(prhs(8))
opts_pr =  mxGetPr(prhs(9))
	
nsx_pr = mxGetPr(plhs(1))
nsy_pr = mxGetPr(plhs(2))


!    load data into Fortran arrays

 call mxCopyPtrToReal8(s_pr,s,nxyz)
 call mxCopyPtrToReal8(t_pr,t,nxyz)
 call mxCopyPtrToReal8(p_pr,p,nxyz)
 call mxCopyPtrToReal8(g_pr,g,nxyz)
 call mxCopyPtrToReal8(ocean_pr,ocean,nxy)
 call mxCopyPtrToReal8(n_pr,n,nxy)
 call mxCopyPtrToReal8(longs_pr,longs,nxy)
 call mxCopyPtrToReal8(lats_pr,lats,nxy)
 call mxCopyPtrToReal8(opts_pr,opts,2)

	
!    call computational subroutine


 call quick_gradients(s,t,p,g,ocean,n,longs,lats,opts,nx,ny,nz,nsx,nsy)
		
 	

!    load output into MATLAB arrays

 call mxCopyReal8ToPtr(nsx,nsx_pr,nxyz)		
 call mxCopyReal8ToPtr(nsy,nsy_pr,nxyz)	

		

return
end     
