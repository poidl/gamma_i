function [ss,tt] = dj_stabilize0(s,t,p,N2)

%%% STABILIZE:      Stabilize hydrographic data w.r.t. buoyancy frequency
%%%
%%% USAGE:          [ss,tt,pp] = dj_stabilize(s,t,p,N2)
%%%
%%% INPUT:          s       matrix of salinity (each column being a cast)
%%%                 t       matrix of in-situ temperatures 
%%%                 p   	 matrix of pressures
%%%                 N2      vector of N2 lower bounds
%%%
%%%                 NOTE:   missing values denoted by NaN's
%%%
%%% OUTPUT:         ss    	 matrix of adjusted salinities
%%%                 tt      matrix of adjusted in-situ temperatures
%%%
%%% UNITS:          salinity    				psu (IPSS-78)
%%%                 temperature degrees 	deg. C (IPTS-68)
%%%                 pressure    				db
%%%                 N2     					sec^-2
%%%
%%%
%%% AUTHOR:         David Jackett
%%%
%%% CREATED:        May, 2000
%%%
%%% REVISION:       1.1     6/3/97
%%%



%%%
%%%     check # arguments and initialize
%%%

if nargin ~= 4
  error('ERROR in dj_stabilize: invalid input arguments')
end

[nz,nx] = size(s);

ss = nan*zeros(nz,nx); tt = ss; pp = ss;

if nz == 1
  np = ~isnan(s+t+p);
else
  np = sum(~isnan(s+t+p));
end

for ix = 1:nx
  index(1:np(ix),ix) = find(~isnan(s(:,ix)+t(:,ix)+p(:,ix)));
end


%%%
%%%     save appropriate array
%%%

savearray = [];
for ix = 1:nx
  indx = index(1:np(ix),ix);
  savearray = [savearray; ...
           along(ix) alat(ix) np(ix); ...
           s(indx,ix) t(indx,ix) p(indx,ix)];
end

zpath = path;

%%	 a mod from Eric Firing <efiring@iniki.soest.hawaii.edu>
%% 		to get the correct path when it is first or last
%%			in the path variable

p_delims = findstr(zpath,';');            % Path delimiter indices.
p_starts = [1 (p_delims+1)];              % Indices of start of each
										  %						subpath.
p_ends   = [(p_delims-1) length(zpath)];  % Indices of end of each subpath.
zgamma = findstr(zpath,'gamma_n');        % Index inside desired subpath.
i_path = max(find (p_starts <= zgamma));  % Desired index into p_starts.

gamma_path = zpath(p_starts(i_path):p_ends(i_path)); % Desired path.

zpwd = pwd;

command = ['cd ', gamma_path]; eval(command);

save glab_m.in savearray -ascii


%%%
%%%     run external code
%%%

command = ['!glab_m']; eval(command);

load glab_m.out;

delete glab_m.in glab_m.out

command = ['cd ',zpwd]; eval(command);


%%%
%%%     assemble gamma_n labels
%%%

start = 1;
for ix = 1:nx
  n = np(ix);
  indx = index(1:n,ix);
  finish = start+n-1;
  g(indx,ix) = glab_m(start:finish,1);
  dg_lo(indx,ix) = glab_m(start:finish,2);
  dg_hi(indx,ix) = glab_m(start:finish,3);
  start = finish+1;
end


return
