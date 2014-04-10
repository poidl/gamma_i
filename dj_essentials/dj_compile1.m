clear all, close all, clc

if exist('d:/software/msdev/LIB/MATHD.LIB')~=0
	default_library = 'd:/software/msdev/LIB/MATHD.LIB';
else
	dj_disp('****  double precision MATHD.LIB library not found  ****')
end


cmd = ['mex dqprog_imsl.f90 ', default_library]; eval(cmd)


dj_disp('compiled successfully'); disp(' ')

