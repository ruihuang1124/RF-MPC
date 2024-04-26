function [Xd_,Ud_] = fcn_gen_JumpXdUd(p)

addpath('../arclab-quad-sdk/nmpc_controller/scripts/utils/casadi-windows/')
import casadi.*
%% parameters
gait = p.gait;

plan_steps = 40;
plan_time_horizon = 0.8;
dt_val = repmat(plan_time_horizon/(plan_steps),1,plan_steps);
Fz_max = 128;
max_jump_z=0.8;%最高跳跃高度  但是不能保证一定跳跃到
min_damp_z=0.15;%最低限制高度
max_lift_spdz=8.5;%最大离地速度
z_init=0.2;%起始站立高度
world.mu=p.mu;%摩擦系数

q_init_val = [0 0 0 ...
              0 0 z_init]';%初始状态 rpy  xyz
qd_init_val = [0 0 0 0 0 0]';

q_term_ref  = [1*pi*0   -22/57.3*0      pi*0  ...
               0.4+1e-5       0          0.1 ]';%终端位置 rpy  xyz
qd_term_ref = [0 0 0, 0 0 0]';

%Phase_duty=[int16(N/Sect*0.2) int16(N/Sect*0.3) int16(N/Sect*0.49) int16(N/Sect*0.01)];%相序百分比 下蹲 蹬腿  腾空  缓冲
cs_val = [repmat([1 1 1 1]', 1, 15)   repmat([1 1 1 1]', 1, 5)    repmat([0 0 0 0]', 1, 19)  repmat([0 0 0 0]', 1, 1)]';%MPCtable相序  

weight.QX = [10 10 10, 10 10 10, 10 10 10, 10 10 10 ]';%系统状态代价的权重 轨迹
weight.QN = [10 10 10, 50 50 150, 10 10 10, 10 10 10 ]';%终端代价
weight.Qc = [0.001 0.001 0.001]';
weight.Qf = [0.0001 0.0001 0.001]';

Body.m = p.mass;
Body.m = 5+0.25*8;%机器人质量
%机身惯量
Body.Ib = p.J;
Body.length_body=p.L;
Body.width_body=p.W;
Body.hipPos=p.pf34;
endPos=[0,         0,          0,          0;%足端跨关节位置  %初始化足端位置
             0,         0,          0,          0;
             -z_init,   -z_init,   -z_init,    -z_init];
endPos_Body=Body.hipPos+endPos;
Phip_swing=reshape(endPos_Body,[],1);
world.g = 9.8;%重力加速度

mu_inv = 1.0/world.mu;
% 摩擦约束
f_block =[ mu_inv, 0,  -1.0;
          -mu_inv, 0,  -1.0;
           0,  mu_inv, -1.0;
           0, -mu_inv, -1.0;];

kin_box_x = 0.15; % 运动学约束
kin_box_y = 0.15;
kin_box_z = 0.3; % 腿最长

Kin_block =[ 1, 0,  0,-kin_box_x; % 髋关节坐标系
            -1, 0,  0,-kin_box_x;
             0, 1,  0,-kin_box_y;
             0, -1, 0,-kin_box_y;
             0, 0,  1,-min_damp_z; % 腿最短
             0, 0, -1,-kin_box_z];

acc_d = p.acc_d;
vel_d = p.vel_d;
yaw_d = p.yaw_d;

%% 构造线性化后的系统微分方程

Xk=SX.sym('Xk', 12, 1);%cassis 符号类
n_state=size(Xk,1);
Fk=SX.sym('Uk', 12, 1);
n_F=size(Fk,1);
Rk=SX.sym('Rk', 12, 1);
n_r=size(Rk,1);

I3=eye(3);%单位矩阵
Rbody=rotsb(Xk(1:3));
cy = cos(Xk(3));
sy = sin(Xk(3));
cp = cos(Xk(2));
sp = sin(Xk(2));

R_yaw =[cy sy 0;
        -sy cy 0;
        0 0 1];%世界到机身
R_w=[cy/cp,sy/cp,0;
    -sy,cy,0;
    cy*sp/cp,sy*sp/cp,1];
Ig = Rbody*Body.Ib*Rbody';
Ig_inv=Ig\I3;

%单刚体动力学模型
A = [zeros(3) zeros(3) R_yaw zeros(3)  ;
     zeros(3) zeros(3) zeros(3) I3 ;
     zeros(3) zeros(3) zeros(3) zeros(3);
     zeros(3) zeros(3) zeros(3) zeros(3) ;
    ];%状态矩阵

% AA=A;
% AA(1:3,7:9)=R_w;
B=[ zeros(3)           zeros(3)           zeros(3)            zeros(3);
    zeros(3)           zeros(3)           zeros(3)            zeros(3);
    Ig_inv*Skew(Rk(1:3)) Ig_inv*Skew(Rk(4:6)) Ig_inv*Skew(Rk(7:9))  Ig_inv*Skew(Rk(10:12));
    I3/Body.m   I3/Body.m   I3/Body.m    I3/Body.m;];%控制矩阵
g=zeros(12,1);
g(12)=-world.g;%扩展一维度重力加速度
dotX=A*Xk+B*Fk+g;%构造微分动力学的符号方程

f = Function('f',{Xk;Fk;Rk},{dotX},{'input_states','control_inputs','foot_input'},{'dotX'});

%% 构造代价和约束
% 变量定义, 系统变量及reference parameters变量

X = SX.sym('X', n_state, plan_steps+1); % N+1步状态
F = SX.sym('F', n_F, plan_steps); % N步内的控制 力控制
r = SX.sym('r', n_r, plan_steps); % N步内的控制 足端位置控制

RefX = SX.sym('RefX', n_state, N+1); % N步内的控制输出
RefF = SX.sym('RefF', n_F, N); % N步内的控制输出
Refr = SX.sym('Refr', n_r, N); % N步内的控制输出
ContactState=SX.sym('ConState', 4, N);

% 变量定义用于构造约束和代价的中间变量定义

defect_init=RefX(:,1)-X(:,1); % 12*1 eq
defect_state=SX.zeros(12*(N+1),1); % 12(N+1)*1 eq
defect_FootOnGround=SX.zeros(4*(N),1); % 4(N)*1 eq
defect_footStance=SX.zeros(12*(N),1); % (3*4)(N)*1 eq

n_equa_c=12+12*(N+1)+4*(N)+12*(N);

defect_legLimits=SX.zeros(24*(N),1);%(4*6)(N)*1
defect_footforce=SX.zeros(16*(N),1);%(4*4)(N)*1 xy摩擦约束4个
defect_ForceNormal=SX.zeros(N,1);% (N)*1
defect_footswing=SX.zeros(4*(N),1);%4(N)*1

n_inequa_c=24*(N)+16*(N)+N+4*(N);

% 约束和代价计算
obj = 0;
for k = 1:plan_steps       
    Xk=X(:,k);
    Fk=F(:,k);
    rk=r(:,k);
    Pk=repmat(Xk(4:6),4,1)+rk;
    ContactStatek=ContactState(:,k);
    dtk=dt_val(k);
    
    X_err = Xk - RefX(:,k);                                         % 基座状态误差
    pf_err = repmat(Xk(4:6),4,1) + Phip_swing - Pk;                 % 悬空时约束foot位置 不影响跳跃  只影响最终足端状态输出的结果
    U_err = Fk - RefF(:,k);                                         % GRF 误差
    obj = obj + (X_err'*diag(weight.QX)*X_err+...                   % 误差求和
          pf_err'*diag(repmat(weight.Qc,4,1))*pf_err+...
          U_err'*diag(repmat(weight.Qf,4,1))*U_err)*dtk;
    defect_state((k-1)*12+1:(k-1)*12+12)=X(:,k+1)-(Xk+f(Xk,Fk,rk)*dtk); % 欧拉积分
    % 法向力大于0 不等式
    defect_ForceNormal(k)=-dot(Fk,repmat([0;0;1],4,1));
    % 结合法向力大于0，摩擦约束来约束摆动中力为0 和最大力 不等式
    defect_footswing((k-1)*4+1:(k-1)*4+4)=Fk([3,6,9,12])-ContactStatek.*repmat(1000,4,1);
    for leg=1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        % 脚在地上约束 0是此时地面高度等式
        defect_FootOnGround((k-1)*4+leg)=ContactStatek(leg)*Pk(3*(leg-1)+3);
        % 限制腿长 限制范围不等式
        Rbody=rotsb(Xk(1:3));
        endWorld=Rbody*endPos_Body+Xk(4:6); % 全局足端位置
        p_rel = (Pk(xyz_idx) - endWorld(:,leg)); % hip->足端  足端矢量
        defect_legLimits((k-1)*24+(leg-1)*6+1:(k-1)*24+(leg-1)*6+6)= Kin_block*[p_rel;1];%运动学约束
        % 接触中脚不滑动
        if (k < plan_steps)
            Pk1=repmat(X(4:6,k+1),4,1)+r(:,k+1);
            defect_footStance((k-1)*12+(leg-1)*3+1:(k-1)*12+(leg-1)*3+3)=ContactStatek(leg)*(Pk1(xyz_idx)-Pk(xyz_idx));
        end
        % 摩擦约束 不等式
        defect_footforce((k-1)*16+(leg-1)*4+1:(k-1)*16+(leg-1)*4+4)=f_block*Fk(xyz_idx);
    end
end
% 约束生成
g=[defect_init;defect_state;defect_FootOnGround;defect_footStance;...
    defect_legLimits;defect_footforce;defect_ForceNormal;defect_footswing];
display_str=['等式约束数量',num2str(n_equa_c),'   不等式约束数量',num2str(n_inequa_c)];
disp(display_str);
% 终端 cost
X_err = X(:,end)-RefX(:,end);
obj = obj + X_err'*diag(weight.QN)*X_err;


%%	根据前面构造的约束，代价函数，reference parameters变量，构造优化问题nlp_prob以及优化问题求解器solver
OPT_variables = [reshape(X,n_state*(plan_steps+1),1);reshape(F,n_F*plan_steps,1);reshape(r,n_r*plan_steps,1)];
OPT_Param = [reshape(RefX,n_state*(plan_steps+1),1);reshape(RefF,n_F*plan_steps,1);reshape(Refr,n_r*plan_steps,1);reshape(ContactState,4*plan_steps,1)];
nlp_prob =struct('f', obj, 'x',OPT_variables,'p',OPT_Param, 'g',g );
%  优化设置
opts_setting.expand =true;
opts_setting.ipopt.max_iter=1500;
opts_setting.ipopt.print_level=0;
opts_setting.ipopt.acceptable_tol=1e-4;
opts_setting.ipopt.acceptable_obj_change_tol=1e-6;
opts_setting.ipopt.tol=1e-4;
opts_setting.ipopt.nlp_scaling_method='gradient-based';
opts_setting.ipopt.constr_viol_tol=1e-3;
opts_setting.ipopt.fixed_variable_treatment='relax_bounds';
% 构造求解器
solver = casadi.nlpsol('solver', 'ipopt', nlp_prob,opts_setting);%可以在线修改足端位置

%% 传入具体约束数值，参考轨迹数值，最终求解优化问题：
%	约束上下界 ，对于等式约束，上下界均为0，对于不等式约束，下界为-inf，上界为0
args.lbg(1:n_equa_c) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_equa_c) = 0;  % 1e-20   % Equality constraints
args.lbg(n_equa_c+1 : n_equa_c+ n_inequa_c) = -inf; % inequality constraints
args.ubg(n_equa_c+1 : n_equa_c+ n_inequa_c) = 0; % inequality constraints
%  状态上边界
tempub=[Body.m*world.g*world.mu*6; Body.m*world.g*world.mu*6 ;Fmax];
args.ubx=[];
UBx=[pi*3*ones(3,1);10*ones(2,1);1; pi*3*ones(3,1);40*ones(3,1)];%状态上界约束跳的最高高度
UBx(6)=max_jump_z;
UBx(12)=max_lift_spdz;
UBx=repmat(UBx,plan_steps+1,1);
UBf=[tempub;tempub;tempub;tempub];
UBf=repmat(UBf,plan_steps,1);
UBp=repmat([0.4;0.4;inf],4,1);
UBp=repmat(UBp,plan_steps,1);
args.ubx=[args.ubx;UBx;UBf;UBp];
%  状态下边界
templb=[-Body.m*world.g*world.mu*6; -Body.m*world.g*world.mu*6 ;0]; %力状态
args.lbx=[];
LBx=[-pi*3*ones(3,1);-10*ones(2,1);min_damp_z; -pi*3*ones(3,1);-40*ones(3,1)];%状态下界
LBx=repmat(LBx,plan_steps+1,1);
LBf=[templb;templb;templb;templb];
LBf=repmat(LBf,plan_steps,1);
LBp=repmat([-0.4;-0.4;-inf],4,1);
LBp=repmat(LBp,plan_steps,1);
args.lbx=[args.lbx;LBx;LBf;LBp];

% 构建跳跃的参考轨迹
c_ref = diag([1 1 1, 1 -1 1, -1 1 1, -1 -1 1])*repmat([0.15 0.094 -z_init],1,4)'; %初始化足端位置
f_ref = zeros(12,1);
% set parameter values 设定期望运动轨迹
for i = 1:6%对状态线性插值
    Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),plan_steps+1);%决定轨迹末端位置
    Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),plan_steps+1);
end
% Z向抛物线
a=[Xref_val(4,1),Xref_val(4,plan_steps/2),Xref_val(4,plan_steps)];%x
b=[q_init_val(6),q_term_ref(6),q_init_val(6)+0.0];%z
Xref_val(6,:) =interp1(a,b,Xref_val(4,:),'spline'); %高度方向做Spline插值
Uref_val=zeros(24,N);
r_ref=zeros(12,N);
for leg = 1:4
    for xyz = 1:3
        Uref_val(3*(leg-1)+xyz,:)= Xref_val(xyz+3,1:end-1) +c_ref(3*(leg-1)+xyz);%F 
        r_ref(3*(leg-1)+xyz,:)= c_ref(3*(leg-1)+xyz);%
        Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N);%P
    end
end

if(1)%线性插值
    for i = 1:6
        Xref_val(i,:)   = linspace(q_init_val(i),q_term_ref(i),N+1);
        Xref_val(6+i,:) = linspace(qd_init_val(i),qd_term_ref(i),N+1);
    end
    for leg = 1:4
        for xyz = 1:3
            Uref_val(3*(leg-1)+xyz,:)    = Xref_val(xyz,1:end-1) + c_ref(3*(leg-1)+xyz);
            Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,N);
        end
    end
end
F_ref=Uref_val(13:24,:);

args.p=[reshape(Xref_val,n_state*(N+1),1);reshape(F_ref,n_F*N,1);reshape(r_ref,n_r*N,1);reshape(cs_val',4*N,1)];%送入了轨迹约束 相序约束
args.x0=[reshape(Xref_val,n_state*(N+1),1);reshape(F_ref,n_F*N,1);reshape(r_ref,n_r*N,1)];%系统初始状态


%% generate reference trajectory
% X = [pc dpc eta wb]
lent = plan_steps;
Xd = zeros(30,lent);
Ud = zeros(12,lent);
Rground = p.Rground;           % ground slope



for ii = 1:lent
    if gait >= 0        % --- March forward and rotate ---
        %%%%%%%%%% linear motion %%%%%%%%%%%%%
        pc_d = [0;0;p.z0];
        dpc_d = [0;0;0];
        for jj = 1:2
            if t(ii) < (vel_d(jj) / acc_d)
                dpc_d(jj) = acc_d * t(ii);
                pc_d(jj) = 1/2 * acc_d * t(ii)^2;
            else
                dpc_d(jj) = vel_d(jj);
                pc_d(jj) = vel_d(jj) * t(ii) - 1/2 * vel_d(jj) * vel_d(jj)/acc_d;
            end
        end
        %%%%%%%%%% angular motion %%%%%%%%%%%%%
        if isempty(Xt)
            ea_d = [0;0;0];
        else
            ea_d = [0;0;yaw_d];
        end
        vR_d = reshape(expm(hatMap(ea_d)),[9,1]);
        wb_d = [0;0;0];
    end
    pfd = reshape(Rground * p.pf34,[12,1]);
    Xd(:,ii) = [pc_d;dpc_d;vR_d;wb_d;pfd];
    
    %%%% force
    if (gait == -3)
        Ud(:,ii) = U_d;
    else
        sum_inStance = sum(bool_inStance(:,ii));
        if sum_inStance == 0    % four legs in swing
            Ud(:,ii) = zeros(12,1);
        else
            Ud([3,6,9,12],ii) = bool_inStance(:,ii)*(p.mass*p.g/sum_inStance);
        end
    end
end

end

%% 工具函数
function rotxm=rotx(theta)
s=sin(theta);
c=cos(theta);
% rotxm=[1,0,0;
%     0,c,s
%     0,-s c]';
rotxm=[1,0,0;
    0,c,-s
    0,s c];
end

function rotym=roty(theta)
s=sin(theta);
c=cos(theta);
% rotym =[c,0,-s;
%     0,1,0;
%     s,0,c]';
rotym =[c,0,s;
    0,1,0;
    -s,0,c];
end

function rotzm=rotz(theta)
s=sin(theta);
c=cos(theta);

% rotzm=[c,s,0;
%     -s,c,0;
%     0,0,1]';
rotzm=[c,-s,0;
    s,c,0;
    0,0,1];
end
%Rsb
function R=rotsb(theta)%构造旋转矩阵
% R=rotx(theta(1))*roty(theta(2))*rotz(theta(3));
R=rotz(theta(3))*roty(theta(2))*rotx(theta(1));

end

function s=Skew(in)
s = [0 -in(3) in(2);
    in(3) 0 -in(1);
    -in(2) in(1) 0];
end