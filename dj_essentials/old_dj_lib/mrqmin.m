function mrq = mrqmin(x, y, sig, a, lista, funcs, par)

% MRQMIN - Levenberg-Marquardt non-linear least squares regression
%
% MRQ = MRQMIN(X, Y, SIG, A, LISTA, FUNCS, PAR)
% Uses Levenberg-Marquardt method, attempting to reduce the chi squared
% value, of a fit between a set of points (X, Y) with individual
% standard deviations SIG (Y and SIG are vectors, X is a vector or matrix
% of length equal to that of Y), and a non-linear function dependent on
% coefficients in vector A.  The vector LISTA numbers the parameters 
% which are actually to be adjusted - the others are fixed at their
% initial values.  The user must supply the function FUNCS, the function
% to be fitted.  FUNCS has the following calling format:
%
% [y, dyda] = funcs(x, a)
%
% where Y is the value of the function at X, given parameters A, and
% DYDA, is the derivative of Y with respect to each parameter A.
% X is assumed to be a row vector (or matrix for multi-variable
% function) of length N, where N is the number of data points.  A is
% a vector of length M, where M is the total number of parameters.
% Y must be returned as a row vector of length N, DYDA must be
% returned as an (M, N) matrix.
%
% Optionally the function can be called with extra parameters, PAR, for any
% specific requirements within FUNCS.  To use PAR, FUNCS has the form:
%
% [y, dyda] = funcs(x, a, par)
%
% The initial call to MRQMIN is required with the fill argument list
% (PAR is optional).  The user then repeat calls to MRQMIN until the
% required convergence is reached.  To perform the next iteration, MRQMIN
% is called as follows:
%
% MRQ = MRQMIN(MRQ)
% where the output of one interation is passed onto the next.
% MRQ is a structure with the following fields (of use to the user):
%
%  yFit     the fitted values of Y.
%  chisq    the chi squared value of the fit.
%  a        the updated best fit parameters.
%
% Typical usage:
%
%    mrq = mrqmin(x, y, sig, a, lista, 'fgauss');
%    while (Convergence test)    % Loop until convergence
%        mrq = mrqmin(mrq);
%    end
%
%    [cv, sd] = covar(mrq);     % Get covariance matrix and uncertainties
%
% Derived from Numerical Recipes
%
% Author: Lindsay Pender
%         CSIRO Marine Research
%         Lindsay.Pender@marine.csiro.au
%
% See also: MRQCOF, COVAR

error(nargchk(1, 7, nargin));
if nargin < 7
   par = [];
end

if nargin > 1
   % Initialization
   % Create our working structure
   
   mrq = struct('x', x, 'y', y, 'sig', sig, 'a', a, 'yFit', [], ...
      'lista', lista, 'funcs', funcs, 'chisq', [], 'lambda', [], 'alpha', [], ...
      'beta', [], 'par', par);
   
   % Check arguments
   
   [ny, my] = size(mrq.y);							% Check for row vector
   if ny * my ~= my
      if ny * my == ny
         mrq.y = mrq.y';
         my = ny;
         
      else
         error('MRQMIN: Y must be a vector');
      end
   end
   
   [n, m] = size(mrq.sig);							% Check for row vector
   if n * m ~= m
      if n * m == n
         mrq.sig = mrq.sig';
         m = n;
         
      else
         error('MRQMIN: SIG must be a vector');
      end
   end
   
   if m ~= my
      error('MRQMIN: SIG must be the same length as X');
   end
   
   [n, m] = size(mrq.x);							% Check for row vector
   if n * m ~= m
      if n * m == n
         mrq.x = mrq.x';
         m = n;
      end
   end
   
   if m ~= my
      error('MRQMIN: X must be the same length as Y');
   end
   
   [na, ma] = size(mrq.a);							% Check for row vector
   if na * ma ~= ma
      if na * ma == na
         mrq.a = mrq.a';
         ma = na;
         
      else
         error('MRQMIN: A must be a vector');
      end
   end
   
   [n, m] = size(mrq.lista);					% Check for row vector
   if n * m ~= m
      if n * m == n
         mrq.lista = mrq.lista';
         m = n;
         
      else
         error('MRQMIN: LISTA must be a vector');
      end
   end
   
   if m > ma
      error('MRQMIN: LISTA contains more elements than A');
   end
   
   i = find(mrq.lista < 1 | mrq.lista > ma);
   if ~isempty(i)
      error('MRQMIN: invalid LISTA elements');
   end
   
   mrq.funcs = funcs;
   if ~ischar(mrq.funcs)
      error('MRQMIN: FUNCS must be a string');
   end
   
   if isempty(which(mrq.funcs))
      error('MRQMIN: FUNCS not found');
   end
   
   mrq.lambda = 0.001;
   [mrq.yFit, mrq.alpha, mrq.beta, mrq.chisq] = mrqcof(mrq, a);
else
   % Next iteration
   
   mrq = x;
   if ~isstruct(mrq)
      error('MRQMIN: MRQ must be a function');
   end
end

% Alter linearized fitting matrix by augmenting diagonal

covar = mrq.alpha;
for j = 1:length(mrq.lista)
   covar(j, j) = covar(j, j) * (1.0 + mrq.lambda);
end

% Get matrix solution

da = mrq.beta / covar;

% Test the trial fit

atry = mrq.a;
for j = 1:length(mrq.lista)
   atry(mrq.lista(j)) = atry(mrq.lista(j)) + da(j);
end

ochisq = mrq.chisq;
[mrq.yFit, covar, da, mrq.chisq] = mrqcof(mrq, atry);
if (mrq.chisq < ochisq)									% Successs - accept solution
   mrq.lambda = mrq.lambda * 0.1;
   mrq.alpha = covar;
   mrq.beta = da;
   mrq.a = atry;
else															% No success - adjust lambda
   mrq.lambda = mrq.lambda * 10;
end
