gn = g;
gi = g;
gn0 = gn;

inds = find(isfinite(gi+gn));
gi = gi(inds);
gn = gn(inds);

%       on raw data
nn = length(gi);
ng = 0;
if ng~=0
    bin_size = floor(nn/ng);
end

gi_adj = nan*ones(ng,1);
gi_nodes = nans(ng,1);

[gi_srtd, k_srtd] = sort(gi);
gn_srtd = gn(k_srtd);

gi_0 = gi_srtd(1);
gi_1 = gi_srtd(nn);

gn_0 = gn_srtd(1);
gn_1 = gn_srtd(nn);

if ng == 0

%    gg = g(:,j0,i0); indsgg = find(finite(gg)); ngg = length(indsgg)    
%    gi_adj = gg(1:ngg); gi_nodes = gi_adj;        
%    inds_g = find(finite(g)); ginds = nan*ones(size(g)); ginds(inds_g) = 1:length(inds_g);   
%    indsmgs = inds_g(ginds(1:ngg,j0,i0));  
%    [indsmgs,g(indsmgs),g(1:ngg,j0,i0),gi_adj]
        
    inds_g = find(isfinite(g));
    ginds = nan*ones(size(g));
    ginds(inds_g) = 1:length(inds_g);
    
%   maximum and minimum casts
    gfac = 0.9;
    
    [gmax,indg] = nanmax(g(30,:));
    ind_gmax = find(g(inds_g) == gmax);

%   [gmax,ind_gmax] = nanmax(g(inds_g)) 
    
    [k0,j0,i0] = ind2sub(size(g),inds_g(ind_gmax(1)));
    
    location = [longs(i0),lats(j0)]
    
    gg = g(:,j0,i0);
    indsgg = find(isfinite(gg));
    ngg = length(indsgg)
    
    indsmgs = inds_g(ginds(1:ngg,j0,i0));
    
    [indsmgs,g(indsmgs),g(1:ngg,j0,i0)]   
    
    [s(1:ngg,j0,i0),t(1:ngg,j0,i0),p(1:ngg,j0,i0),g(1:ngg,j0,i0)]    
    
    gmid = (gg(1)+gmax)/2
    
    indsgg = find(gg>=(gmax-gfac*gmid));
    ngg = length(indsgg)

    indsmgs = inds_g(ginds(1:ngg,j0,i0));
    
    gi_adj = g(indsmgs);
    gi_nodes = gi_adj;   
    
    [indsmgs,g(indsmgs),g(1:ngg,j0,i0),gi_adj]    
 
    k0 = -1;
    gg = g;
    
    while k0~=1
        if k0~=-1
            gg(k0,j0,i0) = nan;
        end
        [gmin,ind_gmin] = nanmin(gg(inds_g))
        [k0,j0,i0] = ind2sub(size(g),inds_g(ind_gmin(1)));
    end
    
    location = [longs(i0),lats(j0)]
    
    gg = g(:,j0,i0);
    indsgg = find(isfinite(gg));
    ngg = length(indsgg)
    
    indsmgs1 = inds_g(ginds(1:ngg,j0,i0));
    
    [indsmgs1,g(indsmgs1),g(1:ngg,j0,i0)]    
      
    indsgg = find(gg<gmin+gfac*gmid);
    ngg = length(indsgg)   
    
    indsmgs = [indsmgs; inds_g(ginds(1:ngg,j0,i0))]; 
    
    gi_adj = [gi_adj; g(inds_g(ginds(1:ngg,j0,i0)))]
    
    [indsmgs,g(indsmgs),gi_adj]    
      
    return
end
    
for k = 1:ng
    inds = [(1+(k-1)*bin_size):(min(nn,k*bin_size))];
    gi_nodes(k) = median(gi_srtd(inds));
    gi_adj(k) = median(gn(k_srtd(inds)));
end

figure
 subplot(2,2,1)
  plot(gi,gn,'.')
  grid on
  hold on
  plot(gi_nodes,gi_adj,'c')
  hold off

 subplot(2,2,2)
  plot(diff(gi_adj))
  grid on
            
%       on 95% data
gi_nodes = [gi_0; gi_nodes; gi_1];
gi_adj = [gn_0; gi_adj; gn_1];
resids = gn - interp1(gi_nodes,gi_adj,gi);
sig = std(resids)

[I95] = find(abs(resids)<=2*sig);
percent_of_data = 100*length(I95)/nn          
gi = gi(I95);
gn = gn(I95);

nn = length(gi);
ng = 0;
bin_size = floor(nn/ng)

gi_nodes = nans(ng,1);
gi_adj = nans(ng,1);

[gi_srtd, k_srtd] = sort(gi);

indsmgs = nans(ng,1);

for k = 1:ng
    inds = [(1+(k-1)*bin_size):(min(nn,k*bin_size))];
    gi_nodes(k) = median(gi_srtd(inds));
    gi_adj(k) = median(gn(k_srtd(inds)));
    [dg_min,indss] = nanmin(abs(gn0(:)-gi_adj(k)));
    indsmgs(k) = indss(1);
end

gi_nodes = [gi_0; gi_nodes; gi_1];
gi_adj = [gn_0; gi_adj; gn_1];

inds0 = find(gn0 == gn_0);
inds0 = inds0(1);

inds1 = find(gn0 == gn_1);
inds1 = inds1(1);

indsmgs = [inds0; indsmgs; inds1];

subplot(2,2,3)
	plot(gi,gn,'.')
    grid on
	hold on
    plot(gi_nodes,gi_adj,'c')
    hold off
            
subplot(2,2,4)
	plot(diff(gi_adj))
    grid on
            
%gi_nodes            
g = interp1(gi_nodes,gi_adj,g);    

