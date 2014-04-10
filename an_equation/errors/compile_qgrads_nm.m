gamman_library =  'd:/mySoftware/gamma_n/fortran/double_precision/gamma_dp.lib';

eos_library = 'd:/eos/eos05/code/fortran/Debug/eoslib.lib';

cmd = ['mex quick_gradients.f90 ns_gradients.f90 sigl_gradients.f90 ', gamman_library, ' ', eos_library];

eval(cmd)

dj_disp('      ... compile_qgrads_nm compiled successfully')

