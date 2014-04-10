function coeffs = generate_coefficients(npoly)


%%              load component coefficients

cd data/R1_nan
    cmd = ['copyfile(''gamma_p', int2str(npoly), '.dat'', ''gamma_p.dat'')']
    eval(cmd), load gamma_p.dat
    c1 = gamma_p;
cd ../R2_nas
    cmd = ['copyfile(''gamma_p', int2str(npoly), '.dat'', ''gamma_p.dat'')']
    eval(cmd), load gamma_p.dat
    c2 = gamma_p;
cd ../R3_saw
    cmd = ['copyfile(''gamma_p', int2str(npoly), '.dat'', ''gamma_p.dat'')']
    eval(cmd), load gamma_p.dat
    c3 = gamma_p;
cd ../R4_sai
    cmd = ['copyfile(''gamma_p', int2str(npoly), '.dat'', ''gamma_p.dat'')']
    eval(cmd), load gamma_p.dat
    c4 = gamma_p;
cd ../..


coeffs = [c1; c2; c3; c4];

           
return