function dj_msave(file,z)

%%%    Usage:    dj_msave(file,z)
%%%
%%%    Description:    save an ascii file of z suitable for MMA
%%%
%%%    Input:          file	- filename (text)
%%%                    z		- array of data to be saved
%%%
%%%    Output:         a file on disk
%%%
%%%    Author:         David Jackett
%%%
%%%    Date:           20/11/96
%%%


cmd = ['exist(''',file,''')']; ifile = eval(cmd);

if ifile == 2, cmd = ['delete(''',file,''')']; eval(cmd); end

cmd = ['save ''',file,''' z -ascii'];

eval(cmd)


return
