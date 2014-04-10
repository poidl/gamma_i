!	parameter(nx=360,ny=161,nz=45)

	parameter(nx=90,ny=43,nz=45)			! Gouretski 4 degree data

	real longs(nx),lats(ny),p(nz)
	real s(nz,ny,nx),t(nz,ny,nx),eta(nz,ny,nx),g(nz,ny,nx)

	real NSx(nz,ny,nx),NSy(nz,ny,nx)
	real SIGLx(nz,ny,nx),SIGLy(nz,ny,nx)
	real SIG0x(nz,ny,nx),SIG0y(nz,ny,nx)
	real GPOLYx(nz,ny,nx),GPOLYy(nz,ny,nx)

	character filename*50, data_set*2


	data_set = 'g4'


	NSq = 1; SIGLq = 1; SIG0q = 0; GPOLYq = 0;


 
    elapsed_time = TIMEF()


!
!						read data
!

	filename = 'data/fdata_'//data_set//'.dat'

	print *, filename

	open(10,file=filename,status='old')


	read(10,*) anx
	read(10,*) any
	read(10,*) anz

	nx1 = anx; ny1 = any; nz1 = anz

	if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
		print *, 'dimensions wrong!'
		stop
	else
		print *, 'dimensions:',nx,ny,nz
	end if
	

	read(10,*) (longs(i),i=1,nx)
	read(10,*) (lats(j),j=1,ny)
	read(10,*) (p(k),k=1,nz)
	read(10,*) (((s(k,j,i),k=1,nz),j=1,ny),i=1,nx)
	read(10,*) (((t(k,j,i),k=1,nz),j=1,ny),i=1,nx)
	read(10,*) (((g(k,j,i),k=1,nz),j=1,ny),i=1,nx)
		
	close(10)

    elapsed_time = TIMEF();
	print *, 'read data: ', elapsed_time,' seconds'


!						
!						compute gradients
!

	if(NSq.eq.1) then
		call NS_gradients(s,t,p,g,longs,lats,nx,ny,nz,NSx,NSy)
	  	open(11,file='data/NS_'//data_set//'.dat',status='unknown')
		write(11,'(e16.8)') NSx,NSy
		close(11)
	end if


	if(SIGLq.eq.1) then
		do j = 1,ny
		do i = 1,nx
		do k = 1,nz
			if(t(k,j,i).gt.-99.0d0) then
				eta(k,j,i) = eta_from_t(s(k,j,i),t(k,j,i),p(k))
			end if
		end do
		end do
		end do
		call SIGL_gradients(s,eta,p,longs,lats,nx,ny,nz,SIGLx,SIGLy)
	  	open(11,file='data/SIGL_'//data_set//'.dat',status='unknown')
		write(11,'(e16.8)') SIGLx,SIGLy
		close(11)
	end if


	if(SIG0q.eq.1) then
		call SIG0_gradients(s,t,p,longs,lats,nx,ny,nz,SIG0x,SIG0y)
	  	open(11,file='data/SIG0_'//data_set//'.dat',status='unknown')
		write(11,'(e16.8)') SIG0x,SIG0y
		close(11)
	end if


	if(GPOLYq.eq.1) then
		call GPOLY_gradients(s,t,p,longs,lats,nx,ny,nz,GPOLYx,GPOLYy)
	  	open(11,file='GPOLYg.dat',status='unknown')
		write(11,'(e16.8)') GPOLYx,GPOLYy
		close(11)
	end if


    elapsed_time = TIMEF();	print *, 'done: ', elapsed_time



	stop
	end