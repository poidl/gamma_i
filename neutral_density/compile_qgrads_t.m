
%gamma_library =  'd:/mySoftware/gamma_n/fortran/double_precision/gamma_dp.lib';

z = matlabroot; 
cmd = ['gamma_library =  ''', z(1),':/mySoftware/gamma_n/fortran/double_precision/gamma_dp.lib'';'];
eval(cmd) 

cmd = ['eos_library =  ''', z(1),':/eos/eos05/code/fortran/Debug/eoslib.lib'';'];
eval(cmd)

%eos_library = 'd:/eos/eos05/code/fortran/Debug/eoslib.lib';

%gfunc_library = './gpoly/gpoly.lib';

sigp_library = './sigp/sigma_p.lib';

cmd = ['mex -v quick_gradients.f90 ns_gradients.f90 sigl_gradients.f90 ', ...
                 eos_library, ' ', gamma_library, ' ', sigp_library];

eval(cmd)

%dj_disp('      compile_qgrads_t compiled successfully ...')

%dj_toc