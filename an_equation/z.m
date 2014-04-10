r0 = 0:0.01:1;  s0 = 0:0.01:1;

[r,s] = meshgrid(r0,s0);

w1 = s; w2 = (1-r).*(1-s); w3 = r.*(1-s);

subplot(2,2,1)
    fpcolor(r0,s0,w1), colorbar
    hold on, contour(r0,s0,w1,20,'k'), hold off
    
subplot(2,2,3)
    fpcolor(r0,s0,w2), colorbar
    hold on, contour(r0,s0,w2,20,'k'), hold off
    
subplot(2,2,2)
    fpcolor(r0,s0,w3), colorbar
    hold on, contour(r0,s0,w3,20,'k'), hold off

wsum = w1+w2+w3;

subplot(2,2,4)
    fpcolor(r0,s0,wsum), colorbar
    hold on, contour(r0,s0,wsum,20,'k'), hold off
    
figure(gcf)
