clear all; close all; clc; dj_tic

pr0 = 0; level = 14;

load s_th

inds = find(pns==-99); sns(inds) = nan; tns(inds) = nan; 

thns = sw_ptmp(sns,tns,pns,pr0);


figure(1)

plot(sns(:),thns(:),'ro'); grid on; hold on

plot(ss,tt,'o')


%		now add the original data

load fig3_stp

ss = reshape(s(level,:),nx,ny); ss = ss';
tt = reshape(t(level,:),nx,ny); tt = tt';
pp = reshape(p(level,:),nx,ny); pp = pp';

inds = find(tt==0); ss(inds) = nan; tt(inds) = nan; pp(inds) = nan;

tthh = sw_ptmp(ss,tt,pp,pr0);

plot(ss(:),tthh(:),'co')

dj_pltmp(longs,lats,tt)


