%dj_tic

gamma_library =  'd:/mySoftware/gamma_n/fortran/double_precision/gamma_dp.lib';

eos_library = 'd:/eos/eos05/code/fortran/Debug/eoslib.lib';

cmd = ['mex quick_depth_ntp.f90 ', gamma_library, ' ', eos_library];

eval(cmd)

dj_disp('      ... compile_qdntp_t compiled successfully')

%dj_toc