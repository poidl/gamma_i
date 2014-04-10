function [N2_i] = gsw_interp_Nsquared(N2,p_mid,p_i)
% gsw_interp_Nsquared                 linear interpolation to p_i on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating 
% variable p_mid. This function finds the values of N2 at p_i on this cast.
%
% VERSION NUMBER: 3.0 (8th April, 2011) 
%
% This fuction was adapted from Matlab's interp1q.
%==========================================================================

p_mid = p_mid(:);

[min_p,Imin_p] = min(p_mid);
[Ishallow] = find(p_i <= min_p);        % Set equal to the shallowest bottle.
if ~isempty(Ishallow)
    N2_i(Ishallow) = N2(Imin_p);
end

[max_p,Imax_p] = max(p_mid);
[Ideep] = find(p_i >= max_p);            % Set equal to the deepest bottle.
if ~isempty(Ideep)
    N2_i(Ideep) = N2(Imax_p);
end

[I] = find(p_i >= min_p & p_i <= max_p);

xi = p_i(I);

x = p_mid;

siz = size(xi);
if ~isscalar(xi)
   [xxi, k] = sort(xi);
   [dum, j] = sort([x;xxi]);
   r(j) = 1:length(j);
   r = r(length(x)+1:end) - (1:length(xxi));
   r(k) = r;
   r(xi==x(end)) = length(x)-1;
   ind = find((r>0) & (r<length(x)));
   ind = ind(:);
   N2_ri = NaN(length(xxi),size(N2,2),superiorfloat(x,N2,xi));
   rind = r(ind);
   xrind = x(rind);
   u = (xi(ind)-xrind)./(x(rind+1)-xrind);
   N2rind = N2(rind,:);
   if exist('bsxfun','builtin') == 5
       N2_ri(ind,:) = N2rind + bsxfun(@times,N2(rind+1,:)-N2rind,u);
   else
       N2_ri(ind,:) = N2rind + (N2(rind+1,:)-n2rind).*u;
   end
else
   % Special scalar xi case
   r = find(x <= xi,1,'last');
   r(xi==x(end)) = length(x)-1;
   if isempty(r) || r<=0 || r>=length(x)
      N2_ri = NaN(1,size(N2,2),superiorfloat(x,N2,xi));
   else
      u = (xi-x(r))./(x(r+1)-x(r));
      N2r = N2(r,:);
      if exist('bsxfun','builtin') == 5
          N2_ri = N2r + bsxfun(@times,N2(r+1,:)-N2r,u);
      else
          N2_ri = N2r + (N2(r+1,:)-N2r).*u;
      end
   end
end

if min(size(N2_ri)) == 1 && numel(xi) > 1
   N2_ri = reshape(N2_ri,siz);
end

N2_i(I) = N2_ri;

end
