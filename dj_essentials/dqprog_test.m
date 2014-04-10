clear all, close all, clc, dj_tic

ncon = 2; nvars = 5; neq = 2;

A = [1 1 1 1 1; 0 0 1 -2 -2]; b = [5 -3]';

g = [-2 0 0 0 0]';

H = [2 0 0 0 0; 0 2 -2 0 0 ; 0 -2 2 0 0 ; 0 0 0 2 -2; 0 0 0 -2 2];

x = dqprog_imsl(A,b,g,H,neq)
