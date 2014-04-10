	parameter(nx=90,ny=43,nz=45)

	implicit real*8(a-h,o-z)

	real*8 longs(nx),lats(ny),p(nz)
	real*8 ocean(ny,nx)
	real*8 s(nz,ny,nx),t(nz,ny,nx),g(nz,ny,nx)
	real*8 nsx(nz,ny,nx),nsy(nz,ny,nx) 
	 
	 
	 
	open(10,file='ocean_dump.dat',status='old')
	read(10,*) nx1,ny1,nz1
	read(10,'(e20.12)') longs,lats,p,s,t,g,ocean
	close(10)

	print *; print *, 'read data'
	 
	call sigp_gradients(s,t,p,longs,lats,nx,ny,nz,nsx,nsy)


	stop
	end