function dj_pause(nsecs)

%%%    dj_pause		modified pause function
%%%
%%%	 Usage:			dj_pause(nsecs)
%%%
%%%	 Input:			nsecs, 0 or # seconds to pause
%%%
%%%    Author:       David Jackett
%%%
%%%    Date:         24/4/96
%%%


if nsecs == 0
  disp('  ok ?   '); pause
else
  pause(nsecs)
end

return