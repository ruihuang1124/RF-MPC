clear all;close all;clc
addpath fcns fcns_MPC
addpath('../qpSWIFT/matlab/')

%% --- parameters ---
% ---- gait ----
% 0-trot; 1-bound; 2-pacing 3-gallop; 4-trot run; 5-crawl; [-6]-complex jump
gait = 4;
p = get_params(gait);
[Xd_,Ut_] = fcn_gen_JumpXdUd(p);
