clear all; close all; clc; dj_tic

x = 0:1:10; y = x; z = x;

[X,Y,Z] = meshgrid(x,y,z);

f = X.*Y.*Y.*Z;

[fx,fy,fz] = gradient(f);

