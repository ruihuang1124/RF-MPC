clear all;close all;clc
addpath fcns fcns_MPC

p.predHorizon = 2;
p.simTimeStep = 1/200;
p.Tmpc = 4/100;             % MPC prediction step time
p.Umax = 50;
p.decayRate = 1;
p.freq = 30;
p.Rground = eye(3);
p.Qf = diag([1e5 2e5 3e5 5e2 1e3 150 1e3 1e4 800 40 40 10]);

p.mass = 5.5;
p.J = diag([0.026,0.112,0.075]);
p.g = 9.81;
p.mu = 1;       % friction coefficient
p.z0 = 0.2;     % nominal COM height
p.pf34 = [[0.15;0.094;0],[0.15;-0.094;0],[-0.15;0.094;0],[-0.15;-0.094;0]];

p.L = 0.301;    % body length
p.W = 0.088;    % body width
p.d = 0.05;     % ABAD offset
p.h = 0.05;     % body height
p.l1 = 0.14;    % link1 length
p.l2 = 0.14;    % link2 length

p.Kp_sw = 300;  % Kp for swing phase

p.body_color    = [42 80 183]/255;
p.leg_color     = [7 179 128]/255;
p.ground_color  = [195 232 243]/255;

p.prejump_time = 0.20;
p.plan_time_horizon = 0.8;
p.plan_steps = p.plan_time_horizon / p.simTimeStep;

p.predHorizon = 2;
p.simTimeStep = 1/100;
p.Tmpc = 0.02;
p.Tst = 0.3;
p.Tsw = 0.15;
p.R = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
p.Q = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
p.Qf = p.Q;
% min. 0.5 * x' * H *x + g' * x
% s.t. Aineq *x <= bineq
%      Aeq * x <= beq
% X = [pc dpc vR wb pf]': [30,1]
% q = [pc dpc eta wb]: [12 1]
% lb/ub - [4,n_hor]


A = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
B = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
d = [0;1;2;3;4;5;6;7;8;9;10;11];
Ut = [0;1;2;3;4;5;6;7;8;9;10;11];

mu = 1;
n_hor = 2;
Umax = 180;
decayRate = 1;

R = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
Q = diag([0 1 2 3 4 5 6 7 8 9 10 11]);
Qf = Q;
[Qx,Qv,Qeta,Qw] = deal(Q(1:3,1:3),Q(4:6,4:6),Q(7:9,7:9),Q(10:12,10:12));
[Qxf,Qvf,Qetaf,Qwf] = deal(Qf(1:3,1:3),Qf(4:6,4:6),Qf(7:9,7:9),Qf(10:12,10:12));

nX = 12;
nU = 12;
Xt = [];
for ii=0:29
    Xt = [Xt;ii];
end
for jj=1:n_hor
    Xd(:,jj)=Xt + jj-1;
    Ud(:,jj)=Ut + jj-1;
end
Xd;
Ud;

%%%%%%% A,B,d matrices for linear dynamics %%%%%%%%%%%
% [A,B,d] = fcn_get_ABD_eta(Xt,Ut,p);

%% Decompose
Rt = reshape(Xt(7:15,1),[3,3]);
qt = [Xt(1:6);[0;0;0];Xt(16:18)];

% lb <= Fz <= ub
Fzd = Ud([3 6 9 12],:);
lb = -1 * Fzd;
ub = 2 * Fzd;

%% Matrices for QP
H = zeros((nX + nU) * n_hor);
g = zeros(size(H,1),1);
Aeq = zeros(nX * n_hor,(nX+nU) * n_hor);
beq = zeros(size(Aeq,1),1);
Aineq_unit = [1 0 -mu;-1 0 -mu;0 1 -mu;0 -1 -mu;0 0 1; 0 0 -1];
nAineq_unit = size(Aineq_unit,1);
Aineq = zeros(4*nAineq_unit*n_hor,(nX+nU)*n_hor);
bineq = zeros(size(Aineq,1),1);
for i_hor = 1:n_hor
    xd = Xd(1:3,i_hor);
    vd = Xd(4:6,i_hor);
    Rd = reshape(Xd(7:15,i_hor),[3,3]);
    wd = Xd(16:18,i_hor);
    
    %% Objective function
    idx_u = (i_hor-1) * (nX + nU) + (1:nU);
    idx_x = (i_hor-1) * (nX + nU) + nU + (1:nX);
    a = Rd' * Rt;
    a1 = logm(Rd' * Rt);
    a2 = veeMap(logm(Rd' * Rt));

    if i_hor == n_hor
        H(idx_x,idx_x) = Qf * decayRate^(i_hor-1);
        g(idx_x) = [-Qxf * xd;
                    -Qvf * vd;
                     Qetaf * veeMap(logm(Rd' * Rt));
                    -Qwf * wd] * decayRate^(i_hor-1);
        b = Qetaf * veeMap(logm(Rd' * Rt));
    else
        H(idx_x,idx_x) = Q * decayRate^(i_hor-1);
        g(idx_x) = [-Qx * xd;
                    -Qv * vd;
                     Qeta * veeMap(logm(Rd' * Rt));
                    -Qw * wd] * decayRate^(i_hor-1);
        b = Qeta * veeMap(logm(Rd' * Rt));
    end
    H(idx_u,idx_u) = R * decayRate^(i_hor-1);
    g(idx_u) = R' * (Ut - Ud(:,i_hor)) * decayRate^(i_hor-1);

                
    %% Equality constraints
    if i_hor == 1
        Aeq(1:nX,1:(nU+nX)) = [-B,eye(nX)];
        beq(1:nX) = A * qt + d;
    else
        Aeq((i_hor-1)*nX+(1:nX),(i_hor-2)*(nX+nU)+nU+(1:(2*nX+nU)))= [-A -B eye(nX)];
        % a = (i_hor-1)*nX+(1:nX)
        % b = (i_hor-2)*(nX+nU)+nU+(1:(2*nX+nU))
        beq((i_hor-1)*nX+(1:nX)) = d;
    end

    %% Inequality constraints
    Fi = zeros(4*nAineq_unit,12);
    hi = zeros(size(Fi,1),1);
    for i_leg = 1:4
        idx_F = (i_leg-1)*nAineq_unit + (1:nAineq_unit);
        idx_u = (i_leg-1)*3 + (1:3);
        Fi(idx_F,idx_u) = Aineq_unit;
        hi(idx_F) = [mu*Ut(idx_u(3))-Ut(idx_u(1));
                     mu*Ut(idx_u(3))+Ut(idx_u(1));
                     mu*Ut(idx_u(3))-Ut(idx_u(2));
                     mu*Ut(idx_u(3))+Ut(idx_u(2));
                     ub(i_leg,i_hor)-Ut(idx_u(3))+Ud(idx_u(3),i_hor);
                     -lb(i_leg,i_hor)+Ut(idx_u(3))-Ud(idx_u(3),i_hor)];
    end
    idx_A = (i_hor-1) * 4*nAineq_unit + (1:4*nAineq_unit);
    idx_z = (i_hor-1) * (nX+nU) + (1:nU);
    Aineq(idx_A,idx_z) = Fi;
    Fi;
    bineq(idx_A) = hi;
end

[Cx_x,Cx_v,Cv_v,Cv_u,Cv_c] = eta_co_xv(Ut,p.Tmpc,p.mass,p.g);
% Cx_x, Cx_v, Cv_v, Cv_u, Cv_c

xop = reshape(Xt(1:3),[3,1]);
vop = reshape(Xt(4:6),[3,1]);
Rop = reshape(Xt(7:15),[3,3]);
wop = reshape(Xt(16:18),[3,1]);
pf34 = reshape(Xt(19:30),[3,4]);
p.Tmpc

[CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,p.Tmpc);
% CE_eta, CE_w, CE_c

[Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,Ut,p.Tmpc,p.J,pf34);
% Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c

[AA,BB,dd] = fcn_get_ABD_eta(Xt,Ut,p);
% AA,BB,dd

[HH,gg,Aineqq,bineqq,Aeqq,beqq] = fcn_get_QP_form_eta(Xt,Ut,Xd,Ud,p);
HH,gg,Aineqq,bineqq,Aeqq,beqq

function [A,B,D] = fcn_get_ABD_eta(Xt,Ut,p)
% linear dynamics for rotation
% evolution variable is eta

%% parameters
dt = p.Tmpc;

%% unpack
xop = reshape(Xt(1:3),[3,1]);
vop = reshape(Xt(4:6),[3,1]);
Rop = reshape(Xt(7:15),[3,3]);
wop = reshape(Xt(16:18),[3,1]);
pf34 = reshape(Xt(19:30),[3,4]);

%% constants for linear matrices
% [x,v,eta,w,constant]
[Cx_x,Cx_v,Cv_v,Cv_u,Cv_c] = eta_co_xv(Ut,dt,p.mass,p.g);
Cx_x,Cx_v,Cv_v,Cv_u,Cv_c
[CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt);
CE_eta, CE_w, CE_c
[Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,Ut,dt,p.J,pf34);
Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c


%% Assemble matrices
    A = [Cx_x, Cx_v, zeros(3,6);
         zeros(3), Cv_v, zeros(3,6);
         zeros(3,6),CE_eta,CE_w;
         Cw_x,zeros(3),Cw_eta,Cw_w];
    B = [zeros(3,12);
         Cv_u;
         zeros(3,12);
         Cw_u];
    D = [zeros(3,1);
         Cv_c;
         CE_c;
         Cw_c];
    A,B,D

end

function [Cx_x,Cx_v,Cv_v,Cv_u,Cv_c] = eta_co_xv(fop,dt,mass,g)

Cx_x = eye(3);
Cx_v = eye(3) * dt;

Cv_v = eye(3);
Cv_u = dt/mass * [eye(3),eye(3),eye(3),eye(3)];
Cv_c = Cv_u * fop + [0;0;-g] * dt;

end

function [CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N

%% debugged code
invN = pinv(N)

C_eta = kron(eye(3),Rop*hatMap(wop))*N + kron(eye(3),Rop)*fcn_get_D(wop);
C_w = kron(eye(3),Rop) * N;
C_c = vec(Rop*hatMap(wop)) - kron(eye(3),Rop)*N*wop;

CE_eta = eye(3) + invN * dt * kron(eye(3),Rop') * C_eta;
CE_w = invN * dt * kron(eye(3),Rop') * C_w;
CE_c = invN * dt * kron(eye(3),Rop') * C_c;

end

function [Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,fop,dt,J,pf)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;
r1 = pf(:,1) - xop;
r2 = pf(:,2) - xop;
r3 = pf(:,3) - xop;
r4 = pf(:,4) - xop;
Mop = [hatMap(r1) hatMap(r2) hatMap(r3) hatMap(r4)] * fop;

temp_J_w = hatMap(J*wop) - hatMap(wop) * J;
sum_fop = [eye(3),eye(3),eye(3),eye(3)] * fop;

Cx = Rop' * hatMap(sum_fop);
Ceta = fcn_get_F(Rop'*Mop) * N - temp_J_w * hatMap(wop);
Cw = temp_J_w;
Cu = Rop' * [hatMap(r1),hatMap(r2),hatMap(r3),hatMap(r4)];
Cc = -hatMap(wop)*J*wop + Rop'*Mop - temp_J_w * wop - Cx*xop;

Cw_x = dt*(J\Cx);
Cw_eta = dt*(J\Ceta);
Cw_w = dt*(J\Cw) + eye(3);
Cw_u = dt*(J\Cu);
Cw_c = dt*(J\Cc);

end

function F = fcn_get_F(k)

F = [k', zeros(1,3),zeros(1,3);...
     zeros(1,3),k',zeros(1,3);...
     zeros(1,3),zeros(1,3),k'];
end

function N = fcn_get_N

N = [0 0 0;...
     0 0 1;...
     0 -1 0;...
     0 0 -1;...
     0 0 0;...
     1 0 0;...
     0 1 0;...
     -1 0 0;...
     0 0 0];
end

function D = fcn_get_D(in)

d = in(1);
e = in(2);
f = in(3);
D = [0 0 0;
   e -d 0;
   f 0 -d;
   -e d 0;
   0 0 0;
   0 f -e;
   -f 0 d;
   0 -f e;
   0 0 0];
end

function [H,g,Aineq,bineq,Aeq,beq] = fcn_get_QP_form_eta(Xt,Ut,Xd,Ud,p)
% min. 0.5 * x' * H *x + g' * x
% s.t. Aineq *x <= bineq
%      Aeq * x <= beq
% X = [pc dpc vR wb pf]': [30,1]
% q = [pc dpc eta wb]: [12 1]
% lb/ub - [4,n_hor]

%% parameters
mu = p.mu;
n_hor = p.predHorizon;
Umax = p.Umax;
decayRate = p.decayRate;

R = p.R;
Q = p.Q;
Qf = p.Qf;
[Qx,Qv,Qeta,Qw] = deal(Q(1:3,1:3),Q(4:6,4:6),Q(7:9,7:9),Q(10:12,10:12));
[Qxf,Qvf,Qetaf,Qwf] = deal(Qf(1:3,1:3),Qf(4:6,4:6),Qf(7:9,7:9),Qf(10:12,10:12));

nX = 12;
nU = 12;

%%%%%%% A,B,d matrices for linear dynamics %%%%%%%%%%%
[A,B,d] = fcn_get_ABD_eta(Xt,Ut,p);

%% Decompose
Rt = reshape(Xt(7:15,1),[3,3]);
qt = [Xt(1:6);[0;0;0];Xt(16:18)];

% lb <= Fz <= ub
Fzd = Ud([3 6 9 12],:);
lb = -1 * Fzd;
ub = 2 * Fzd;

%% Matrices for QP
H = zeros((nX + nU) * n_hor);
g = zeros(size(H,1),1);
Aeq = zeros(nX * n_hor,(nX+nU) * n_hor);
beq = zeros(size(Aeq,1),1);

Aineq_unit = [1 0 -mu;-1 0 -mu;0 1 -mu;0 -1 -mu;0 0 1; 0 0 -1];

nAineq_unit = size(Aineq_unit,1);
Aineq = zeros(4*nAineq_unit*n_hor,(nX+nU)*n_hor);
bineq = zeros(size(Aineq,1),1);
for i_hor = 1:n_hor
    xd = Xd(1:3,i_hor);
    vd = Xd(4:6,i_hor);
    Rd = reshape(Xd(7:15,i_hor),[3,3]);
    wd = Xd(16:18,i_hor);
    
    %% Objective function
    idx_u = (i_hor-1) * (nX + nU) + (1:nU);
    idx_x = (i_hor-1) * (nX + nU) + nU + (1:nX);
    if i_hor == n_hor
        H(idx_x,idx_x) = Qf * decayRate^(i_hor-1);
        g(idx_x) = [-Qxf * xd;
                    -Qvf * vd;
                     Qetaf * veeMap(logm(Rd' * Rt));
                    -Qwf * wd] * decayRate^(i_hor-1);
    else
        H(idx_x,idx_x) = Q * decayRate^(i_hor-1);
        g(idx_x) = [-Qx * xd;
                    -Qv * vd;
                     Qeta * veeMap(logm(Rd' * Rt));
                    -Qw * wd] * decayRate^(i_hor-1);
    end
    H(idx_u,idx_u) = R * decayRate^(i_hor-1);
    g(idx_u) = R' * (Ut - Ud(:,i_hor)) * decayRate^(i_hor-1);

                
    %% Equality constraints
    if i_hor == 1
        Aeq(1:nX,1:(nU+nX)) = [-B,eye(nX)];
        beq(1:nX) = A * qt + d;
    else
        Aeq((i_hor-1)*nX+(1:nX),(i_hor-2)*(nX+nU)+nU+(1:(2*nX+nU)))= [-A -B eye(nX)];
        beq((i_hor-1)*nX+(1:nX)) = d;
    end

    %% Inequality constraints
    Fi = zeros(4*nAineq_unit,12);
    hi = zeros(size(Fi,1),1);
    for i_leg = 1:4
        idx_F = (i_leg-1)*nAineq_unit + (1:nAineq_unit);
        idx_u = (i_leg-1)*3 + (1:3);
        Fi(idx_F,idx_u) = Aineq_unit;
        % if p.gait == -2
        %     hi(idx_F) = [Umax-Ut(idx_u(1));Umax+Ut(idx_u(1));...
        %                  Umax-Ut(idx_u(2));Umax+Ut(idx_u(2));...
        %                  Umax-Ut(idx_u(3));Umax+Ut(idx_u(3))];
        % else
        %     hi(idx_F) = [mu*Ut(idx_u(3))-Ut(idx_u(1));
        %                  mu*Ut(idx_u(3))+Ut(idx_u(1));
        %                  mu*Ut(idx_u(3))-Ut(idx_u(2));
        %                  mu*Ut(idx_u(3))+Ut(idx_u(2));
        %                  ub(i_leg,i_hor)-Ut(idx_u(3))+Ud(idx_u(3),i_hor);
        %                 -lb(i_leg,i_hor)+Ut(idx_u(3))-Ud(idx_u(3),i_hor)];
        % end
        hi(idx_F) = [mu*Ut(idx_u(3))-Ut(idx_u(1));
                     mu*Ut(idx_u(3))+Ut(idx_u(1));
                     mu*Ut(idx_u(3))-Ut(idx_u(2));
                     mu*Ut(idx_u(3))+Ut(idx_u(2));
                     ub(i_leg,i_hor)-Ut(idx_u(3))+Ud(idx_u(3),i_hor);
                     -lb(i_leg,i_hor)+Ut(idx_u(3))-Ud(idx_u(3),i_hor)];
    end
    idx_A = (i_hor-1) * 4*nAineq_unit + (1:4*nAineq_unit);
    idx_z = (i_hor-1) * (nX+nU) + (1:nU);
    Aineq(idx_A,idx_z) = Fi;
    bineq(idx_A) = hi;
end



end
