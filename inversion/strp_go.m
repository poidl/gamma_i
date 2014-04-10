clear all, close all, clc

load andreas

disp('converting to pp ...')
    tic, Ap = matlab2pp(A);
         bp = matlab2pp(b); toc
         
disp('setting up ...')         
    tic, AtAp = Ap'*Ap;
         Atbp = Ap'*bp; toc
         
disp('inverting ...')
    tic, gp = AtAp\Atbp; toc
    
disp('writing output ...')
    tic, gg = pp2matlab(gp); 
         save andreas_gg gg
         toc

