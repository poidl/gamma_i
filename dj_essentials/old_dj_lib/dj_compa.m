function dj_compa(array1,array2)

%%%    dj_compa		compare two 2-dimensional arrays
%%%
%%%    Usage:    dj_compa(array1,array2)
%%%
%%%    Description:    compare two arrays
%%%
%%%    Input: 			  array1 - first array
%%%						  array2 - second array
%%%
%%%    Output:         prints the rms difference of the two arrays,
%%%					     and when this is nonzero, a histogram of the
%%%				        differences
%%%
%%%    Author:         David Jackett
%%%
%%%    Date:           12/4/96


inds = find(finite(array1+array2));

darray = array1(inds)-array2(inds);

dmin = min(darray); dmax = max(darray);

rms = sqrt(sum(darray.*darray))/length(darray);

disp(['  {', num2str(dmin), ' min, ' num2str(rms), ' rms, ', ...
			num2str(dmax), ' max} difference'])

if rms ~= 0
  hist(darray,50)
end

clear inds darray dmin dmax rms

return
