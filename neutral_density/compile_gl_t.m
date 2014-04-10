
z = matlabroot; 

cmd = ['gamma_library =  ''', z(1),':/mySoftware/gamma_n/fortran/double_precision/gamma_dp.lib'';'];
eval(cmd)

cmd = ['eos_library =  ''', z(1),':/eos/eos05/code/fortran/Debug/eoslib.lib'';'];
eval(cmd)

cmd = ['mex quick_section_glabel.f90 ', gamma_library, ' ', eos_library];

eval(cmd)

cmd = ['mex quick_global_glabel.f90 ', gamma_library, ' ', eos_library];

eval(cmd)

% dj_disp('      ...   compiled successfully ...')
