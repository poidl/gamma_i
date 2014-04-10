global s t ct p g longs lats

longss = longs(1,:); latss = lats(:,1);

[x,y] = ginput(1);

[xmin,i] = min(abs(x-longss)); [ymin,j] = min(abs(y-latss));

%j = 102; i = 99;

location = [i,longss(i),j,latss(j)]

% inds = find(finite(g(:,j,i)));


ss = s(:,j,i); ctt = ct(:,j,i); pp = p(:,j,i); gg = g(:,j,i);


% [gg1,dggl,dggh] = gamma_n(ss,tt,pp,longs(i),lats(j));

% NN2 = bfrq(ss,tt,pp,'ct');
% 
% ctt = ct_from_t(ss,tt,pp);

data = [ss,ctt,pp,gg]

% 
% 
% figure, subplot(2,2,1), plot(ss,ctt,'.-'), grid on
%             title(['(',int2str(round(longss(i))),'\circE,',int2str(round(latss(j))),'\circN)'])
%             xlabel('\itS'), ylabel('\Theta')
%         subplot(2,2,3),  plot(gg,pp), grid on, set(gca,'ydir','reverse')
%             xlabel('\gamma'), ylabel('\itp')
%         subplot(2,2,4),  plot([NN2;nan],pp), grid on, set(gca,'ydir','reverse')
%             xlabel('N^2'), ylabel('\itp')
% 
% 
