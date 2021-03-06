function [siglx,sigly] = ntp_gradients(s,t,p,g,ocean,n,longs,lats)

[nz,ny,nx] = size(s); plevel = 200;


siglx = nan*ones(nz,ny,nx); sigly = siglx;

s_east = nan*ones(nz,1); t_east = s_east;
s_west = nan*ones(nz,1); t_west = s_west;
s_north = nan*ones(nz,1); t_north = s_north;
s_south = nan*ones(nz,1); t_south = s_south;

sns_e = nan*ones(nz,1); tns_e = sns_e; pns_e = sns_e;
sns_w = nan*ones(nz,1); tns_w = sns_w; pns_w = sns_w;
sns_n = nan*ones(nz,1); tns_n = sns_n; pns_n = sns_n;
sns_s = nan*ones(nz,1); tns_s = sns_s; pns_s = sns_s;

for j0 = 2:ny-1

for i0 = 2:nx-1

  denx = 111.2d3*cos(pi*lats(j0,i0)/180)*(longs(j0,i0+1)-longs(j0,i0-1));
  deny = 111.2d3*(lats(j0+1,i0)-lats(j0-1,i0));
  if ocean(j0,i0)>=1 & ocean(j0,i0)<=7 & ...
       s(1,j0,i0)>=0 & s(1,j0,i0+1)>=0 & ...	
 	     s(1,j0+1,i0)>=0 & s(1,j0,i0-1)>=0 & ...	
   	      s(1,j0-1,i0)>=0

    n_east = 0;
    for k1 = 1:nz
      if s(k1,j0,i0+1)>=0
        n_east = n_east+1; 
        s_east(n_east) = s(k1,j0,i0+1);
        t_east(n_east) = t(k1,j0,i0+1);
        p_east(n_east) = p(k1,j0,i0+1);
      end
    end

    n_west = 0;
    for k1 = 1:nz
      if s(k1,j0,i0-1)>=0
        n_west = n_west+1; 
        s_west(n_west) = s(k1,j0,i0-1);
        t_west(n_west) = t(k1,j0,i0-1);
        p_west(n_west) = p(k1,j0,i0-1);
      end
    end

    n_north = 0;
    for k1 = 1:nz
      if s(k1,j0+1,i0)>=0
        n_north = n_north+1; 
        s_north(n_north) = s(k1,j0+1,i0);
        t_north(n_north) = t(k1,j0+1,i0);
        p_north(n_north) = p(k1,j0+1,i0);
      end
    end

    n_south = 0;
    for k1 = 1:nz
      if s(k1,j0-1,i0)>=0
        n_south = n_south+1; 
        s_south(n_south) = s(k1,j0-1,i0);
        t_south(n_south) = t(k1,j0-1,i0);
        p_south(n_south) = p(k1,j0-1,i0);
      end
    end


%		now the work
    
    for k0 = 1:nz
      if p(k0,j0,i0)>=plevel & s(k0,j0,i0)>=0
        [sns_e(k0),tns_e(k0),pns_e(k0)] = ...
            depth_ntp(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0),s_east(1:n_east),t_east(1:n_east),p_east(1:n_east));
        [sns_w(k0),tns_w(k0),pns_w(k0)] = ...
            depth_ntp(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0),s_west(1:n_west),t_west(1:n_west),p_west(1:n_west));
        [sns_n(k0),tns_n(k0),pns_n(k0)] = ...
            depth_ntp(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0),s_north(1:n_north),t_north(1:n_north),p_north(1:n_north));
        [sns_s(k0),tns_s(k0),pns_s(k0)] = ...
            depth_ntp(s(k0,j0,i0),t(k0,j0,i0),p(k0,j0,i0),s_south(1:n_south),t_south(1:n_south),p_south(1:n_south));
      else
        pns_e(k0) = -99;
        pns_w(k0) = -99;
        pns_n(k0) = -99;
        pns_s(k0) = -99;
      end
    end

%          the x-gradient

    for k0 = 1:nz
      if p(k0,j0,i0)>=plevel & pns_e(k0)>=0 & pns_w(k0)>=0
        siglx(k0,j0,i0) =  -(pns_e(k0)-pns_w(k0))/denx;
      end
    end

%          the y-gradient

    for k0 = 1:nz
      if p(k0,j0,i0)>=plevel & pns_n(k0)>=0 & pns_s(k0)>=0
        sigly(k0,j0,i0) =  -(pns_n(k0)-pns_s(k0))/deny;
      end
    end

  end

end
j0, dj_toc
end

return
