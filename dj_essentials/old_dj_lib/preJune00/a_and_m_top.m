
%%					arctic top
 	
	s = reshape(s1(1,:),ny,nx);
	th = reshape(th1(1,:),ny,nx);
  p = reshape(p1(1,:),ny,nx);
	pr0 = p_arctic*ones(ny,nx);
%%        %%%%%%
	t = sw_ptmp(s,th,pr0,p);
	sig0 = sw_pden(s,t,p,pr0)-1000;
 
	inds = find((ocean==7|ocean==9)&sig0>=sig_arctic);
  nzeros = length(inds);
%%                               %%%%%%
	pns(level,inds) = p1(1,1)*ones(nzeros,1);


%%					mediterannean top
 	
	s = reshape(s1(1,:),ny,nx);
	th = reshape(th1(1,:),ny,nx);
  p = reshape(p1(1,:),ny,nx);
	pr0 = p_med*ones(ny,nx);
%%        %%%%%%
	t = sw_ptmp(s,th,pr0,p);
	sig0 = sw_pden(s,t,p,pr0)-1000;
 
	inds = find((ocean==8)&sig0>=sig_med);
  nzeros = length(inds);
%%                               %%%%%%
	pns(level,inds) = p1(1,1)*ones(nzeros,1);


