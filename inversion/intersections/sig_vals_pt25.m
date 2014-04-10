function [sig1,sig2] =  sig_vals_pt25(SP1,pt1,p1,SP2,pt2,p2)

%	Computes the sigma values of two neighbouring bottles w.r.t. the mid pressure
%
%   Usage :         [sig1,sig2] =  sig_vals(SP1,pt1,p1,SP2,pt2,p2)
%
%	Input :			SP1, SP2         bottle salinities
%					pt1, pt2         bottle Conservative Temperatures
%					p1, p2           bottle pressures
%
%	Output :		sig1, sig2       bottle potential density values
%
%	Units :			SP      g/kg (EOS-80)
%					pt      degrees C (IPTS-90)
%					pressure      dbar
%                   density       kg m-3

pmid = 0.5*(p1 + p2);
sig1 = rho_from_pt25(SP1,pt1,pmid) - 1000;
sig2 = rho_from_pt25(SP2,pt2,pmid) - 1000;

end