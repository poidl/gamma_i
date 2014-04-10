dj_tic


%%              load rfunc_7_9 coefficients

cd data/R1a
    copyfile('rfunc_7_9.dat','rfunc.dat'), load rfunc.dat
    c1 = rfunc

cd ../R2a
    copyfile('rfunc_7_9.dat','rfunc.dat'), load rfunc.dat
    c2 = rfunc
    
cd ../..


%%              

[nz,ny,nx] = size(g);

wt1 = nan*ones(size(n)); wt2 = wt1; wt3 = wt1; wt4 = wt1;

inds_g = find(finite(g(1,:))); wt1(inds_g) = 0; wt2(inds_g) = 0; wt3(inds_g) = 0; wt4(inds_g) = 0;

grf = nan*ones(size(g));

lat_lo = 32; lat_hi = 40;

for j = 1:ny
    for i = 1:nx
        
        if ocean(j,i)==5 & finite(g(1,j,i))
            
            if lats(j,i)>=lat_hi
                wt1(j,i) = 1;
            elseif lats(j,i)<=lat_lo
                wt2(j,i) = 1;
            else
                wt1(j,i) = (lats(j,i)-lat_lo)/(lat_hi-lat_lo);  wt2(j,i) = 1-wt1(j,i);
            end  
            
            c = wt1(j,i)*c1 + wt2(j,i)*c2;
            
            indsz = find(finite(g(:,j,i))); %[indsz,s(indsz,j,i),ct(indsz,j,i)]
            
            gg = rfunc_7_9(s(indsz,j,i),ct(indsz,j,i),c);
            
            grf(indsz,j,i) = gg;

        end
    end
end

g = grf-1000;

% longss = nanmean(longs); latss = nanmean(lats'); 
% 
% figure(1), subplot(2,2,1), dj_pltmp(longss,latss,wt1), title('wt1')
% 
%            subplot(2,2,2), dj_pltmp(longss,latss,wt2), title('wt2')
% 
%             subplot(2,2,3), dj_pltmp(longss,latss,ocean)
% 
%            subplot(2,2,4), dj_pltmp(longss,latss,wt4), title('wt4')
% 
% 
% colormap(flipud(jet))



%%      and plot

smin = nanmin(s(:)); smax = nanmax(s(:)); sby = (smax-smin)/100;

ctmin = nanmin(ct(:)); ctmax = nanmax(ct(:)); ctby = (ctmax-ctmin)/100;

s0 = smin:sby:smax; ct0 = ctmin:ctby:ctmax; ns = length(s0); nct = length(ct0);

[ss,ctt] = meshgrid(s0,ct0); ss = ss(:); ctt = ctt(:);

g0 = rfunc_7_9(ss,ctt,c); 

g0 = reshape(g0,nct,ns)-1000;

figure(5)
    fpcolor(s0,ct0,g0), colorbar
    inds = find(finite(s(:))); ss = s(inds); ctt = ct(inds);
    hold on, plot(ss,ctt,'w.')
    contour(s0,ct0,g0,20,'k'), hold off
    
dj_toc