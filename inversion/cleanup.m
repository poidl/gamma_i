global s t ct p g longs lats ocean n
inds = find(g<=0);
length(inds)
[nz,ny,nx] = size(g);
%[indsz,indsy,indsx] = ind2sub(size(g),inds);
%figure(1), hold on, subplot(2,2,3)
%plot(longs(indsx),lats(indsy),'mx'), hold off
g(inds) = nan; s(inds) = nan; t(inds) = nan; ct(inds) = nan; p(inds) = nan;