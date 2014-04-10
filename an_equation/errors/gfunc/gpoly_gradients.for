	subroutine  gpoly_gradients(s,t,p,longs,lats,nx,ny,nz,
     &				gpolyx,gpolyy)

	parameter(nlev=45)


	real longs(nx),lats(ny),p(nz)
	real s(nz,ny,nx),t(nz,ny,nx)

	real gpoly(nlev,189,192)

	real s_east(nlev), t_east(nlev), p_east(nlev), gpoly_east(nlev)
	real s_west(nlev), t_west(nlev), p_west(nlev), gpoly_west(nlev)
	real s_north(nlev),t_north(nlev),p_north(nlev),gpoly_north(nlev)
	real s_south(nlev),t_south(nlev),p_south(nlev),gpoly_south(nlev)

	real pns_e(nlev),pns_w(nlev),pns_n(nlev),pns_s(nlev)

	real gpolyx(nz,ny,nx),gpolyy(nz,ny,nx)

	data pi/3.14159265/, pr0/0.0/, plevel/200.d0/



c
c						initialize variables
c

	gpolyx = -99.0; gpolyy = -99.0; ng = 1;

	do j0 = 1,ny
	do i0 = 1,nx
	do k0 = 1,nz
		if(s(k0,j0,i0).gt.0.0) then
			th0 = theta(s(k0,j0,i0),t(k0,j0,i0),p(k0),pr0)
			gpoly(k0,j0,i0) = dj_gamma(s(k0,j0,i0),th0)
		end if
	end do
	end do
	end do



c
c						main loop
c		

	do j0 = 2,ny-1
	deny = 111.2e3*(lats(j0+1)-lats(j0-1))
	do i0 = 2,nx-1

		denx = 111.2e3*cos(pi*lats(j0)/180)*(longs(i0+1)-longs(i0-1))

		if(s(1,j0,i0).ne.-99.9.and.s(1,j0,i0+1).ne.-99.9.and.
     &		s(1,j0+1,i0).ne.-99.9.and.s(1,j0,i0-1).ne.-99.9.and.
	&									s(1,j0-1,i0).ne.-99.9) then

			n_east = 0
			do k1 = 1,nz
				if(s(k1,j0,i0+1).gt.0.0) then
					n_east = n_east+1 
					s_east(n_east) = s(k1,j0,i0+1);
					t_east(n_east) = t(k1,j0,i0+1);
					p_east(n_east) = p(k1);
					gpoly_east(n_east) = gpoly(k1,j0,i0+1);
				end if
			end do

			n_west = 0
			do k1 = 1,nz
				if(s(k1,j0,i0-1).gt.0.0) then
					n_west = n_west+1 
					s_west(n_west) = s(k1,j0,i0-1);
					t_west(n_west) = t(k1,j0,i0-1);
					p_west(n_west) = p(k1);
					gpoly_west(n_west) = gpoly(k1,j0,i0-1);
				end if
			end do

			n_north = 0
			do k1 = 1,nz
				if(s(k1,j0+1,i0).gt.0.0) then
					n_north = n_north+1 
					s_north(n_north) = s(k1,j0+1,i0);
					t_north(n_north) = t(k1,j0+1,i0);
					p_north(n_north) = p(k1);
					gpoly_north(n_north) = gpoly(k1,j0+1,i0);
				end if
			end do

			n_south = 0
			do k1 = 1,nz
				if(s(k1,j0-1,i0).gt.0.0) then
					n_south = n_south+1 
					s_south(n_south) = s(k1,j0-1,i0);
					t_south(n_south) = t(k1,j0-1,i0);
					p_south(n_south) = p(k1);
					gpoly_south(n_south) = gpoly(k1,j0-1,i0);
				end if
			end do


c															now the work

			do k0 = 1,nz

			  if(p(k0).ge.plevel.and.s(k0,j0,i0).gt.0.0) then

				call gpoly_surfaces(s_east,t_east,p_east,gpoly_east,
     &							n_east,gpoly(k0,j0,i0),pns_e(k0))

				call gpoly_surfaces(s_west,t_west,p_west,gpoly_west,
     &							n_west,gpoly(k0,j0,i0),pns_w(k0))

			call gpoly_surfaces(s_north,t_north,p_north,gpoly_north,
     &							n_north,gpoly(k0,j0,i0),pns_n(k0))

			call gpoly_surfaces(s_south,t_south,p_south,gpoly_south,
     &							n_south,gpoly(k0,j0,i0),pns_s(k0))

			  else
				pns_e(k0) = -99.9
				pns_w(k0) = -99.9
				pns_n(k0) = -99.9
				pns_s(k0) = -99.9
	          end if

		  	end do

c						the x-gradient

			do k0 = 1,nz
				if(p(k0).ge.plevel.and.pns_e(k0).gt.0.and.
     &											pns_w(k0).gt.0.0) then
					gpolyx(k0,j0,i0) =  -(pns_e(k0)-pns_w(k0))/denx
				else
					gpolyx(k0,j0,i0) = -99.9
				end if
			end do

c						the y-gradient

			do k0 = 1,nz
				if(p(k0).ge.plevel.and.pns_n(k0).gt.0.and.
     &											pns_s(k0).gt.0.0) then
					gpolyy(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny
				else
					gpolyy(k0,j0,i0) = -99.9
				end if
			end do

		end if

	end do
	end do


	
	return
	end
