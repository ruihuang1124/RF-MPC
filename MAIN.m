% Representation-Free Model Predictive Control for Dynamic Quadruped Panther
% Author: Yanran Ding
% Last modified: 2020/12/21
% 
% Code accompanying the paper:
% Yanran Ding, Abhishek Pandala, Chuanzheng Li, Young-Ha Shin, Hae-Won Park
% "Representation-Free Model Predictive Control for Dynamic Motions in Quadrupeds"
% Transactions on Robotics
% 
% preprint available at: https://arxiv.org/abs/2012.10002
% video available at: https://www.youtube.com/watch?v=iMacEwQisoQ&t=101s

%% initialization
clear all;close all;clc
addpath fcns fcns_MPC
addpath('../qpSWIFT/matlab/')

%% --- parameters ---
% ---- gait ----
% 0-trot; 1-bound; 2-pacing 3-gallop; 4-trot run; 5-crawl; [-6]-complex jump
gait = 0;
p = get_params(gait);
p.playSpeed = 5;
p.flag_movie = 0;       % 1 - make movie
use_qpSWIFT = 1;        % 0 - quadprog, 1 - qpSWIFT (external)
optimal_jump = 1;

dt_sim = p.simTimeStep;
SimTimeDuration = 2.8;  % [sec]
MAX_ITER = floor(SimTimeDuration/p.simTimeStep);

% desired trajectory
p.acc_d = 1;
p.vel_d = [0.5;0];
p.yaw_d = 0;

%% Model Predictive Control
% --- initial condition ---
% Xt = [pc dpc vR wb pf]': [30,1]
if gait == 1
    [p,Xt,Ut] = fcn_bound_ref_traj(p);
else
    [Xt,Ut] = fcn_gen_XdUd(0,[],[1;1;1;1],p);
    if optimal_jump == 1
        [Xd_,Ud_] = fcn_gen_JumpXdUd(p);
    end
end

% --- logging ---
tstart = 0;
tend = dt_sim;


[tout,Xout,Uout,Xdout,Udout,Uext,FSMout] = deal([]);

% --- simulation ----
h_waitbar = waitbar(0,'Calculating...');
tic
for ii = 1:MAX_ITER
    % --- time vector ---
    t_ = dt_sim * (ii-1) + p.Tmpc * (0:p.predHorizon-1);
    if optimal_jump == 1
        if ii < (p.plan_steps)
            Xd_(:,ii+p.plan_steps) = Xd_(:,p.plan_steps);
            Ud_(:,ii+p.plan_steps) = Ud_(:,p.plan_steps);
            %loading optimized trajectory
            FSM = zeros(4,1);
            Xd = Xd_(:,ii:ii+p.predHorizon - 1);
            Ud = Ud_(:,ii:ii+p.predHorizon - 1);
            Ut_test = Ud_(:,ii);
        else
            t_ = t_ - dt_sim * (p.plan_steps-1);

            % [Xt_test,Ut_test] = fcn_gen_XdUd(0,[],[1;1;1;1],p);
            % 
            % for jj = 1:p.predHorizon
            %     Xd(:,jj) = Xd_(:,p.plan_steps);
            %     Ud(:,jj) = Ud_(:,p.plan_steps);
            % end
            
          % --- FSM ---
            if gait == 1
                [FSM,Xd,Ud,Xt] = fcn_FSM_bound(t_,Xt,p);
            else
                [FSM,Xd,Ud,Xt] = fcn_FSM(t_,Xt,p);
            end
       end
    else
      % --- FSM ---
        if gait == 1
           [FSM,Xd,Ud,Xt] = fcn_FSM_bound(t_,Xt,p);
         else
           [FSM,Xd,Ud,Xt] = fcn_FSM(t_,Xt,p);
         end
    end
    % --- MPC ----
    % form QP
    [H,g,Aineq,bineq,Aeq,beq] = fcn_get_QP_form_eta(Xt,Ut,Xd,Ud,p);

    if ~use_qpSWIFT
        % solve QP using quadprog
        [zval] = quadprog(H,g,Aineq,bineq,Aeq,beq,[],[]);
    else
        % interface with the QP solver qpSWIFT
        [zval,basic_info] = qpSWIFT(sparse(H),g,sparse(Aeq),beq,sparse(Aineq),bineq);
    end
    
    % Ut = Ut + zval(1:12);
    % Ut_test = Ud(:,1);
    if ii >= p.plan_steps
        Ut = Ut + zval(1:12);
        Ut_test = Ut;
    end
    
    % Ut = Ut + zval(1:12);
    % Ut_test = Ut;
    % --- external disturbance ---
    [u_ext,p_ext] = fcn_get_disturbance(tstart,p);
    p.p_ext = p_ext;        % position of external force
    u_ext = 0*u_ext;
    
    % --- simulate ---
    % [t,X] = ode45(@(t,X)dynamics_SRB(t,X,Ut,Xd,0*u_ext,p),[tstart,tend],Xt);
    [t,X] = ode45(@(t,X)dynamics_SRB(t,X,Ut_test,Xd,0*u_ext,p),[tstart,tend],Xt);

    
    % --- update ---
    Xt = X(end,:)';
    tstart = tend;
    tend = tstart + dt_sim;
    
    % --- log ---  
    lent = length(t(2:end));
    tout = [tout;t(2:end)];
    Xout = [Xout;X(2:end,:)];
    % Uout = [Uout;repmat(Ut',[lent,1])];
    Uout = [Uout;repmat(Ut_test',[lent,1])];
    Xdout = [Xdout;repmat(Xd(:,1)',[lent,1])];
    Udout = [Udout;repmat(Ud(:,1)',[lent,1])];
    Uext = [Uext;repmat(u_ext',[lent,1])];
    FSMout = [FSMout;repmat(FSM',[lent,1])];
    
    waitbar(ii/MAX_ITER,h_waitbar,'Calculating...');
end
close(h_waitbar)
fprintf('Calculation Complete!\n')
toc

%% Animation
[t,EA,EAd] = fig_animate(tout,Xout,Uout,Xdout,Udout,Uext,p);





