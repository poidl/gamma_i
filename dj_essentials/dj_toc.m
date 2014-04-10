function dj_toc

%%%    dj_toc            modified toc function
%%%
%%%    Author:           David Jackett
%%%
%%%    Date:             24/11/96


t = toc;

if t < 60
  t_string = ['  {', num2str(t,3), ' secs}'];
elseif t < 3600
  t_string = ['  {', num2str(t/60,3), ' mins}'];
else
  t_string = ['  {', num2str(t/3600), ' hrs}'];
end


disp(t_string)

disp(' ')
