function [yFit, alpha, beta, chisq] = mrqcof(mrq, a)

% MRQCOF - compute components for MRQMIN
%
% Used by MRQMIN to evaluate the latest fitted values, the linearized
% fitting matrix ALPHA, vector BETA and CHISQ.
%
% Derived from Numerical Recipes
%
% Author: Lindsay Pender
%         CSIRO Marine Research

nData = length(mrq.x);
nParam = length(a);
nFit = length(mrq.lista);

alpha = zeros(nFit, nFit);
beta = zeros(1, nFit);
chisq = 0.0;

if isempty(mrq.par)
   [yFit, dyda] = eval([mrq.funcs, '(mrq.x, a)']);
else
   [yFit, dyda] = eval([mrq.funcs, '(mrq.x, a, mrq.par)']);
end

sig2 = 1.0 ./ (mrq.sig .* mrq.sig);
dy = mrq.y - yFit;
for j = 1:nFit
   wt = dyda(mrq.lista(j), :) .* sig2;
	for k = 1:j
      alpha(j, k) = sum(wt .* dyda(mrq.lista(k), :));
   end
   
	beta(j) = sum(dy .* wt);
end
      
chisq = sum(dy .* dy .* sig2);

% Fill in the symmetric side

for j=2:nFit
   alpha(1:j-1, j) = alpha(j, 1:j-1)';
end
