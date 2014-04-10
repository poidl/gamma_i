subroutine dqprog_imsl(ncons,nvars,A,b,g,H,neq,x)
	
implicit none 

integer*4                             :: ncons, nvars, neq, nact
integer*4, dimension(:), allocatable  :: iact

real*8                                :: diag
real*8,    dimension(ncons)           :: b
real*8,    dimension(nvars)           :: g, x
real*8,    dimension(:), allocatable  :: alamda
real*8,    dimension(ncons,nvars)     :: A
real*8,    dimension(nvars,nvars)     :: H

external dqprog



allocate(iact(nvars)); allocate(alamda(nvars))

call dqprog(nvars,ncons,neq,A,ncons,b,g,H,nvars,diag,x,nact,iact,alamda)


return
end


!----------------------------------------------------------------------------!


!     the matlab mexfile


subroutine mexFunction(nlhs,plhs,nrhs,prhs)
	
implicit none 

integer*4                               :: nlhs,nrhs
integer*4                               :: A_pr,b_pr,g_pr,H_pr,neq_pr,x_pr
integer*4                               :: mxGetM,mxGetN,mxGetPr,mxCreateFull
integer*4                               :: nvars,ncons,neq,size1,size2,n1
integer*4, dimension(*)                 :: plhs, prhs

real*8                                  :: eneq
real*8,    dimension(:),   allocatable  :: b,g,x
real*8,    dimension(:,:), allocatable  :: A,H 



!     check number of arguments

if(nrhs.ne.5) then
  call mexErrMsgTxt('5 inputs required.')
elseif(nlhs.ne.1) then
  call mexErrMsgTxt('1 output required.')
endif


!     sizes of the input matrices

ncons = mxGetM(prhs(1))
nvars = mxGetN(prhs(1))

size1 = ncons*nvars
size2 = nvars*nvars
n1 = 1

allocate(b(ncons)); allocate(g(nvars)); allocate(x(nvars))

allocate(A(ncons,nvars)); allocate(H(nvars,nvars))


!     check to see dimensions balance

if (size1.ne.ncons*nvars) then
  open(10,file='here.txt',status='unknown')
  write(10,*) ncons,nvars
  close(10)
  call mexErrMsgTxt('*** dimensions do not balance ... see here file ***')
endif


!     create matrices for the return arguments

plhs(1) = mxCreateFull(nvars,n1,0)


!     set up pointers

A_pr = mxGetPr(prhs(1))
b_pr = mxGetPr(prhs(2))
g_pr = mxGetPr(prhs(3))
H_pr = mxGetPr(prhs(4))
neq_pr = mxGetPr(prhs(5))

x_pr = mxGetPr(plhs(1))


!     load data into Fortran arrays

call mxCopyPtrToReal8(A_pr,A,size1)
call mxCopyPtrToReal8(b_pr,b,ncons)
call mxCopyPtrToReal8(g_pr,g,nvars)
call mxCopyPtrToReal8(H_pr,H,size2)
call mxCopyPtrToReal8(neq_pr,eneq,n1)

neq = eneq

!     call the computational subroutine

call dqprog_imsl(ncons,nvars,A,b,g,H,neq,x)


!     load the output into a Matlab array

call mxCopyReal8ToPtr(x,x_pr,nvars)


return
end