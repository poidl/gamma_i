function result = dj_rms(v)

%%
% function rms = dj_rms(v)
%
% v                  vector
%
% result             rms error of v

result = norm(v)/sqrt(length(v));
