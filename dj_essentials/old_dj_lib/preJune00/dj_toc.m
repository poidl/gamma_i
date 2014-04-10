function dj_toc

%%%    dj_toc            modified toc function
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             24/11/96


t = toc; nf = flops;

if t < 60
  t_string = ['  {', num2str(t,3), ' secs, '];
elseif t < 3600
  t_string = ['  {', num2str(t/60,3), ' mins, '];
else
  t_string = ['  {', num2str(t/3600), ' hrs, '];
end

if nf < 1.e3
  fl_string = [num2str(nf,4), ' flops}'];
elseif nf < 1.e6
  fl_string = [num2str(nf/1.e3,4), ' Kflops}'];
elseif nf < 1.e9
  fl_string = [num2str(nf/1.e6,4), ' Mflops}'];
else
  fl_string = [num2str(nf/1.e9,4), ' Gflops}'];
end

disp([t_string, fl_string])

disp(' ')
