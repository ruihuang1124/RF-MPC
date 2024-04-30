function [Xd_,Ud_] = fcn_gen_JumpXdUd(p)

addpath('../arclab-quad-sdk/nmpc_controller/scripts/utils/casadi-windows/')
import casadi.*
%% parameters
t_plot = (0:p.simTimeStep:p.plan_time_horizon);
plan_steps = p.plan_steps;
plan_time_horizon = p.plan_time_horizon;
dt_steps = repmat(plan_time_horizon/(plan_steps),1,plan_steps);
dt = p.simTimeStep;
max_force_z = 128;
max_jump_z=0.8; 
number_of_legs = 4;
world.mu=p.mu;
world.g = p.g;
min_position_z=0.15;
max_jump_speed_z=3.5;
prejump_steps = round((p.prejump_time / p.plan_time_horizon) * plan_steps); 

position_z_init=0.2; 
v_z_init = 0.0;
position_z_takeoff = position_z_init;
v_z_takeoff = 0.0;
position_z_max = 0.5;
v_z_maxz = 0.0;
position_z_touch_down = position_z_init;
v_z_touch_down = 0.0;
position_z_final = position_z_init;
v_z_final = 0.0;

position_x_init = 0.0;
v_x_init = 0.0;
position_x_takeoff = position_x_init;
v_x_takeoff = 0.0;
position_x_touch_down = 0.2;
v_x_touch_down = 0.0;
position_x_final = position_x_touch_down;
v_x_final = 0.0;
position_x_maxz = (position_x_touch_down - position_x_takeoff)/2;
v_x_maxz = 0.0;

position_y_init = 0.0;
v_y_init = 0.0;
position_y_takeoff = position_y_init;
v_y_takeoff = 0.0;
position_y_touch_down = 0.1;
v_y_touch_down = 0.0;
position_y_final = position_y_touch_down;
v_y_final = 0.0;
position_y_maxz = (position_y_touch_down - position_y_takeoff)/2;
v_y_maxz = 0.0;


q_init_value = [0 0 0 position_x_init position_y_init position_z_init]'; % rpy  xyz
qd_init_value = [0 0 0 0 0 0]';

% 根据最高点z值以及降落点x值计算抛物线，确定takeoff起飞时候和tauchdown降落时候的状态
% 根据抛物线先确定takeoff的value
v_z_takeoff = sqrt(2*world.g * (position_z_max - position_z_takeoff));
t_flight = 2*v_z_takeoff/world.g;
flight_steps = round((t_flight / p.plan_time_horizon) * plan_steps);
v_x_takeoff = (position_x_touch_down - position_x_takeoff) * world.g / (2*v_z_takeoff);
v_y_takeoff = (position_y_touch_down - position_y_takeoff) * world.g / (2*v_z_takeoff);
qd_takeoff_value = [0 0 0, v_x_takeoff, v_y_takeoff, v_z_takeoff]';
q_takeoff_value = [0 0 0, position_x_takeoff, position_y_takeoff, position_z_takeoff]';
% 抛物线最高点value；
v_x_maxz = v_x_takeoff;
v_y_maxz = v_y_takeoff;
q_max_z_value = [0 0 0, position_x_maxz position_y_maxz position_z_max]'; % rpy  xyz
qd_max_z_value = [0 0 0, v_x_maxz v_y_maxz 0]';
% 触地时刻的value，主要是速度，z向有重力加速度，x向则无加速度
v_x_touch_down = v_x_maxz;
v_y_touch_down = v_y_maxz;
v_z_touch_down = -v_z_takeoff;
q_touchdown_value = [0 0 0, position_x_touch_down, position_y_touch_down, position_z_touch_down]';
qd_touchdown_value = [0 0 0, v_x_touch_down, v_y_touch_down, v_z_touch_down]';
% 
% 最后时刻的value
q_final_reference_value  = [0, 0, 0, position_x_final, position_y_final, position_z_init]';
qd_final_reference_value = [0 0 0, 0 0 0]';%最终全部速度为0

touchdown_steps = plan_steps - prejump_steps - flight_steps;
% predefine carefully:
contact_states_value = [repmat([1 1 1 1]', 1, prejump_steps)  ... % 1: init--->takeoff
    repmat([0 0 0 0]', 1, flight_steps) ...  % 2: takeoff--->touchdown.
    repmat([1 1 1 1]', 1, touchdown_steps)]'; % 3: touchdown--->final 

weight.Q_states = [10 10 10, 10 10 10, 10 10 10, 10 10 10 ]'; 
weight.Q_final_state = [10 10 10, 50 50 150, 10 10 10, 10 10 100 ]';
weight.Q_foot_positions = [0.001 0.001 0.001]';
weight.Q_control_input = [0.0001 0.0001 0.001]';

robot.m = p.mass;
robot.Ib = p.J;
robot.length_body=p.L;
robot.width_body=p.W;
robot.hipPos=p.pf34;
com_to_foot_vector_z_robot_frame_init=[0,         0,          0,          0; 
             0,         0,          0,          0;
             -position_z_init,   -position_z_init,   -position_z_init,    -position_z_init];
foot_pose_under_robot_frame_init=robot.hipPos+com_to_foot_vector_z_robot_frame_init;
foot_position_swing_under_robot_frame=reshape(foot_pose_under_robot_frame_init,[],1);
mu_inv = 1.0/world.mu;
frition_cone =[ mu_inv, 0,  -1.0;
          -mu_inv, 0,  -1.0;
           0,  mu_inv, -1.0;
           0, -mu_inv, -1.0;];

kinematic_box_x = 0.15; 
kinematic_box_y = 0.15;
kinematic_box_z = 0.25; 

kinematic_block =[ 1, 0,  0,-kinematic_box_x; % under hip frame
            -1, 0,  0,-kinematic_box_x;
             0, 1,  0,-kinematic_box_y;
             0, -1, 0,-kinematic_box_y;
             0, 0,  1,-min_position_z;
             0, 0, -1,-kinematic_box_z];

%% Constructing linearized system differential equations

state_X_kstep = SX.sym('Xk', 12, 1); % CasADi SX type.
number_of_states = size(state_X_kstep,1);
control_inputs_U_kstep = SX.sym('Uk', 12, 1);
number_of_control_inputs = size(control_inputs_U_kstep,1);
com_to_foot_vector_world_frame_R_kstep = SX.sym('Rk', 12, 1);
number_of_com_to_foot_vector_R_kstep = size(com_to_foot_vector_world_frame_R_kstep,1);

I3 = eye(3);

rotation_matrix_world_to_robot = rotsb(state_X_kstep(1:3));
cy = cos(state_X_kstep(3));
sy = sin(state_X_kstep(3));
% cp = cos(state_X_kstep(2));
% sp = sin(state_X_kstep(2));

rotation_matrix_world_to_robot_yaw =[cy sy 0;
        -sy cy 0;
        0 0 1]; % world to robot
Ig = rotation_matrix_world_to_robot*robot.Ib*rotation_matrix_world_to_robot';
Ig_inv=Ig\I3;

% SRBD
A = [zeros(3) zeros(3) rotation_matrix_world_to_robot_yaw zeros(3)  ;
     zeros(3) zeros(3) zeros(3) I3 ;
     zeros(3) zeros(3) zeros(3) zeros(3);
     zeros(3) zeros(3) zeros(3) zeros(3) ;
    ];

% AA=A;
% AA(1:3,7:9)=R_w;

B = [ zeros(3)           zeros(3)           zeros(3)            zeros(3);
    zeros(3)           zeros(3)           zeros(3)            zeros(3);
    Ig_inv*Skew(com_to_foot_vector_world_frame_R_kstep(1:3)) Ig_inv*Skew(com_to_foot_vector_world_frame_R_kstep(4:6)) Ig_inv*Skew(com_to_foot_vector_world_frame_R_kstep(7:9))  Ig_inv*Skew(com_to_foot_vector_world_frame_R_kstep(10:12));
    I3/robot.m   I3/robot.m   I3/robot.m    I3/robot.m;];

g_vector = zeros(12,1);
g_vector(12) = -world.g;
dot_state_X = A*state_X_kstep+B*control_inputs_U_kstep+g_vector; % Constructing symbolic equations for differential dynamics

system_dynamic = Function('system_dynamic',{state_X_kstep; control_inputs_U_kstep; com_to_foot_vector_world_frame_R_kstep},{dot_state_X},{'input_states','control_inputs','foot_input'},{'dotX'});

%% Construction cost function and constraints
% Variable definitions, system variables and reference parameters variables

states_X = SX.sym('X', number_of_states, plan_steps+1); % plan_steps+1
control_inputs_U = SX.sym('U', number_of_control_inputs, plan_steps); 
com_to_foot_vector_world_frame_R = SX.sym('r', number_of_com_to_foot_vector_R_kstep, plan_steps);

reference_states_X = SX.sym('RefX', number_of_states, plan_steps+1); 
reference_control_inputs_U = SX.sym('RefU', number_of_control_inputs, plan_steps); 
reference_com_to_foot_vector_world_frame_R = SX.sym('RefR', number_of_com_to_foot_vector_R_kstep, plan_steps); 
contact_states = SX.sym('ConState', number_of_legs, plan_steps);

% Variable definitions are used to construct intermediate variable definitions for constraints and costs.
% defect_** for constraints and error_** for costs calculation.

defect_init = reference_states_X(:,1)-states_X(:,1); % 12*1 eq
defect_state = SX.zeros(12*(plan_steps+1),1); % 12(N+1)*1 eq
defect_FootOnGround = SX.zeros(4*(plan_steps),1); % 4(N)*1 eq
defect_footStance = SX.zeros(12*(plan_steps),1); % (3*4)(N)*1 eq

number_of_equation_constraints = 12+12*(plan_steps+1)+4*(plan_steps)+12*(plan_steps);

defect_legLimits = SX.zeros(24*(plan_steps),1); %(4*6)(N)*1
defect_footforce = SX.zeros(16*(plan_steps),1); %(4*4)(N)*1: xy friction constraints with number of 4.
defect_ForceNormal = SX.zeros(plan_steps,1); % (N)*1
defect_footswing = SX.zeros(4*(plan_steps),1);% 4(N)*1

number_of_inequation_constraints = 24*(plan_steps)+16*(plan_steps)+plan_steps+4*(plan_steps);

% Constraints and cost calculations
cost_object = 0;
for k = 1:plan_steps       
    state_X_current = states_X(:,k);
    control_inputs_U_current = control_inputs_U(:,k);
    com_to_foot_vector_world_frame_R_current = com_to_foot_vector_world_frame_R(:,k);
    foot_position_under_world_frame_Pk_current = repmat(state_X_current(4:6),4,1) + com_to_foot_vector_world_frame_R_current;
    contact_states_current = contact_states(:,k);
    dt_current_step = dt_steps(k);
    
    error_states_X = state_X_current - reference_states_X(:,k);                         % 
    error_foot_position = repmat(state_X_current(4:6),4,1) + foot_position_swing_under_robot_frame - foot_position_under_world_frame_Pk_current;  % 
    error_control_inputs_U = control_inputs_U_current - reference_control_inputs_U(:,k);  % GRF
    cost_object = cost_object + (error_states_X'*diag(weight.Q_states)*error_states_X+...           % sum
          error_foot_position'*diag(repmat(weight.Q_foot_positions,4,1))*error_foot_position+...
          error_control_inputs_U'*diag(repmat(weight.Q_control_input,4,1))*error_control_inputs_U)*dt_current_step;
    defect_state((k-1)*12+1:(k-1)*12+12)=states_X(:,k+1)-(state_X_current+system_dynamic(state_X_current,control_inputs_U_current,com_to_foot_vector_world_frame_R_current)*dt_current_step); % 欧拉积分
    % normal force of f_z > 0
    defect_ForceNormal(k)=-dot(control_inputs_U_current,repmat([0;0;1],4,1));
    defect_footswing((k-1)*4+1:(k-1)*4+4) = control_inputs_U_current([3,6,9,12])-contact_states_current.*repmat(1000,4,1);
    for leg=1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        % foot on the ground constraints
        defect_FootOnGround((k-1)*4+leg) = contact_states_current(leg)*foot_position_under_world_frame_Pk_current(3*(leg-1)+3);
        % leg length limits
        rotation_matrix_world_to_robot = rotsb(state_X_current(1:3));
        desired_foot_position_world_frame = rotation_matrix_world_to_robot*foot_pose_under_robot_frame_init+state_X_current(4:6); % 全局足端位置
        p_rel = (foot_position_under_world_frame_Pk_current(xyz_idx) - desired_foot_position_world_frame(:,leg)); % hip->足端  足端矢量
        defect_legLimits((k-1)*24+(leg-1)*6+1:(k-1)*24+(leg-1)*6+6)= kinematic_block*[p_rel;1]; % 运动学约束
        % no slippery when contact
        if (k < plan_steps)
            Pk1=repmat(states_X(4:6,k+1),4,1)+com_to_foot_vector_world_frame_R(:,k+1);
            defect_footStance((k-1)*12+(leg-1)*3+1:(k-1)*12+(leg-1)*3+3)=contact_states_current(leg)*(Pk1(xyz_idx)-foot_position_under_world_frame_Pk_current(xyz_idx));
        end
        % friction constraints
        defect_footforce((k-1)*16+(leg-1)*4+1:(k-1)*16+(leg-1)*4+4)=frition_cone*control_inputs_U_current(xyz_idx);
    end
end
% 约束生成
g=[defect_init;defect_state;defect_FootOnGround;defect_footStance;...
    defect_legLimits;defect_footforce;defect_ForceNormal;defect_footswing];
display_str=['number of equal constraints: ',num2str(number_of_equation_constraints),' number of inequal constraints: ',num2str(number_of_inequation_constraints)];
disp(display_str);
% 终端 cost
error_states_X = states_X(:,end)-reference_states_X(:,end);
cost_object = cost_object + error_states_X'*diag(weight.Q_final_state)*error_states_X;


%%	According to the previously constructed constraints, 
% cost function, reference parameters variables,
% construct the optimization problem nlp_prob and the optimization problem solver solver
OPT_variables = [reshape(states_X,number_of_states*(plan_steps+1),1);reshape(control_inputs_U,number_of_control_inputs*plan_steps,1);reshape(com_to_foot_vector_world_frame_R,number_of_com_to_foot_vector_R_kstep*plan_steps,1)];
OPT_Param = [reshape(reference_states_X,number_of_states*(plan_steps+1),1);reshape(reference_control_inputs_U,number_of_control_inputs*plan_steps,1);reshape(reference_com_to_foot_vector_world_frame_R,number_of_com_to_foot_vector_R_kstep*plan_steps,1);reshape(contact_states,4*plan_steps,1)];
nlp_prob =struct('f', cost_object, 'x',OPT_variables,'p',OPT_Param, 'g',g );
%  setting for optimal solver
opts_setting.expand =true;
opts_setting.ipopt.max_iter=1500;
opts_setting.ipopt.print_level=0;
opts_setting.ipopt.acceptable_tol=1e-4;
opts_setting.ipopt.acceptable_obj_change_tol=1e-6;
opts_setting.ipopt.tol=1e-4;
opts_setting.ipopt.nlp_scaling_method='gradient-based';
opts_setting.ipopt.constr_viol_tol=1e-3;
opts_setting.ipopt.fixed_variable_treatment='relax_bounds';
% contructing optimal solver
solver = casadi.nlpsol('solver', 'ipopt', nlp_prob,opts_setting);

%% Pass in the specific constraint values, refer to the trajectory values, and finally solve the optimization problem:
%	Constraint upper and lower bounds. 
%   For equality constraints, the upper and lower bounds are both 0. 
%   For inequality constraints, the lower bound is -inf and the upper bound is 0.
args.lbg(1:number_of_equation_constraints) = 0;  % -1e-20  % Equality constraints
args.ubg(1:number_of_equation_constraints) = 0;  % 1e-20   % Equality constraints
args.lbg(number_of_equation_constraints+1 : number_of_equation_constraints+ number_of_inequation_constraints) = -inf; % inequality constraints
args.ubg(number_of_equation_constraints+1 : number_of_equation_constraints+ number_of_inequation_constraints) = 0; % inequality constraints
% upper bound of state
tempub=[robot.m*world.g*world.mu*6; robot.m*world.g*world.mu*6 ;max_force_z];
args.ubx=[];
UBx=[pi*3*ones(3,1);10*ones(2,1);1; pi*3*ones(3,1);40*ones(3,1)]; % The upper bound of the state constrains the maximum height of the jump
UBx(6)=max_jump_z;
UBx(12)=max_jump_speed_z;
UBx=repmat(UBx,plan_steps+1,1);
UBf=[tempub;tempub;tempub;tempub];
UBf=repmat(UBf,plan_steps,1);
UBp=repmat([0.4;0.4;inf],4,1);
UBp=repmat(UBp,plan_steps,1);
args.ubx=[args.ubx;UBx;UBf;UBp];
% lower bound of state
templb=[-robot.m*world.g*world.mu*6; -robot.m*world.g*world.mu*6 ;0]; % 
args.lbx=[];
LBx=[-pi*3*ones(3,1);-10*ones(2,1);min_position_z; -pi*3*ones(3,1);-40*ones(3,1)]; %
LBx=repmat(LBx,plan_steps+1,1);
LBf=[templb;templb;templb;templb];
LBf=repmat(LBf,plan_steps,1);
LBp=repmat([-0.4;-0.4;-inf],4,1);
LBp=repmat(LBp,plan_steps,1);
args.lbx=[args.lbx;LBx;LBf;LBp];

% Construct a reference trajectory for a jump
% 起始站姿作为参考的 在世界坐标系下表示的 从com--feet的向量
com_to_foot_vector_under_world_frame_ref = diag([1 1 1, 1 -1 1, -1 1 1, -1 -1 1])*repmat([0.15 0.094 -position_z_init],1,4)'; % 初始化足端位置

f_ref = zeros(12,1);
% set parameter values 设定期望运动轨迹


for i = 1:6 % 对状态进行插值
    %首先是prejump的阶段，这时完全是触地的，主要生成takeoff 时候的期望速度
    Xref_val(i,1:prejump_steps) = linspace(q_init_value(i),q_takeoff_value(i),prejump_steps); %
    Xref_val(6+i,1:prejump_steps) = linspace(qd_init_value(i),qd_takeoff_value(i),prejump_steps);
    %然后是飞行阶段，这时完全是腾空的，全由takeoff的状态和重力加速度g决定,先线性插值，后面对z方向单独进行处理
    Xref_val(i,prejump_steps+1:prejump_steps + flight_steps) = linspace(q_takeoff_value(i),q_touchdown_value(i),flight_steps);
    Xref_val(6+i,prejump_steps+1:prejump_steps + flight_steps) = linspace(qd_takeoff_value(i),qd_touchdown_value(i),flight_steps);
    %最后是触地阶段，这时完全是触地的，主要前往最终的期望姿态即稳定在速度和加速度均为0的状态；
    Xref_val(i,prejump_steps+flight_steps+1:plan_steps+1) = linspace(q_touchdown_value(i),q_final_reference_value(i),plan_steps + 1 - prejump_steps - flight_steps);
    Xref_val(6+i,prejump_steps+flight_steps+1:plan_steps+1) = linspace(qd_touchdown_value(i),qd_final_reference_value(i),plan_steps + 1 - prejump_steps - flight_steps);
end
% 飞行阶段Z向抛物线,插值速度以及位置
for step = (prejump_steps + 1) : (prejump_steps + flight_steps)
    Xref_val(12,step) = v_z_takeoff - world.g * dt * (step - prejump_steps); % v = v_0 - gt
    Xref_val(6,step) = position_z_takeoff + v_z_takeoff *  dt * (step - prejump_steps) - 0.5 * world.g *  (dt * (step - prejump_steps))^2; % h = h_init + v_0 - 0.5gt^2;
end

Uref_val = zeros(24,plan_steps); %12 feet_position_under_world_frame & 12 force;
r_ref = zeros(12,plan_steps);
for leg = 1:4
    for xyz = 1:3
        Uref_val(3*(leg-1)+xyz,:) = Xref_val(xyz+3,1:end-1) + com_to_foot_vector_under_world_frame_ref(3*(leg-1)+xyz); % position of feet
        r_ref(3*(leg-1)+xyz,:) = com_to_foot_vector_under_world_frame_ref(3*(leg-1)+xyz); % reference com_to_foot_vector_under_world_frame;
        Uref_val(12+3*(leg-1)+xyz,:) = f_ref(xyz).*ones(1,plan_steps); % force of feet
    end
end

F_ref=Uref_val(13:24,:);

args.p=[reshape(Xref_val,number_of_states*(plan_steps+1),1);reshape(F_ref,number_of_control_inputs*plan_steps,1);reshape(r_ref,number_of_com_to_foot_vector_R_kstep*plan_steps,1);reshape(contact_states_value',4*plan_steps,1)];%送入了轨迹约束 相序约束
args.x0=[reshape(Xref_val,number_of_states*(plan_steps+1),1);reshape(F_ref,number_of_control_inputs*plan_steps,1);reshape(r_ref,number_of_com_to_foot_vector_R_kstep*plan_steps,1)]; % 系统初始状态

sol=solver('x0',args.x0,'lbx', args.lbx,'ubx', args.ubx,'lbg', args.lbg,'ubg', args.ubg,'p',args.p); % 调用求解器,输入数据

x_li=sol.x(1:number_of_states*(plan_steps+1));
x_li=reshape(full(x_li),number_of_states,(plan_steps+1)); %

f_li=sol.x(number_of_states*(plan_steps+1)+1 : number_of_states*(plan_steps+1)+number_of_control_inputs*plan_steps);
f_li=reshape(full(f_li),number_of_control_inputs,plan_steps); % force

r_li=sol.x(number_of_states*(plan_steps+1)+number_of_control_inputs*plan_steps+1 : number_of_states*(plan_steps+1)+number_of_control_inputs*plan_steps+number_of_com_to_foot_vector_R_kstep*plan_steps);
r_li=reshape(full(r_li),number_of_control_inputs,plan_steps); % com_to_foot_world_frame
p_li=r_li+repmat(x_li(4:6,1:end-1),4,1); % foot_position_world_frame

figure(11)
plot(t_plot,Xref_val(6,:),'b--',...
    t_plot,x_li(6,:),'b',...
    t_plot,Xref_val(5,:),'r--',...
    t_plot,x_li(5,:),'r',...
    t_plot,Xref_val(4,:),'g--',...
    t_plot,x_li(4,:),'g',...
    'linewidth',1);
xlabel('t(s)') 
ylabel('value(m)')
title('reference & optimized positions "blue for z, green for x, and red for y"')
hold on; grid on;

figure(18)
plot(t_plot(1:end-1),f_li(3,:),'b',...
    t_plot(1:end-1),f_li(1,:),'g',...
    t_plot(1:end-1),f_li(2,:),'r',...
    'linewidth',1);
xlabel('t(s)') 
ylabel('value(N)')
title('optimized control force for each leg "blue for z, green for x, and red for y"')
hold on; grid on;

for i=1:plan_steps
    cube_animate(x_li(:,i),i,p_li(:,i),~contact_states_value(i,:),[0;0;0;0],...
        f_li(:,i),3,[],[],[],[],[],[-20,14],dt_steps,[]);
pause(0.01); %
end

%% generate reference trajectory
Xd_ = zeros(30,plan_steps);
Ud_ = zeros(12,plan_steps);
Rground = p.Rground;           % ground slope
for ii = 1:plan_steps
    pf34 = [[p_li(1,ii);p_li(2,ii);p_li(3,ii)],[p_li(4,ii);p_li(5,ii);p_li(6,ii)],[p_li(7,ii);p_li(8,ii);p_li(9,ii)],[p_li(10,ii);p_li(11,ii);p_li(12,ii)]];
    pc_d = [x_li(4,ii);x_li(5,ii);x_li(6,ii)];
    dpc_d = [x_li(10,ii);x_li(11,ii);x_li(12,ii)];
    wb_d = [x_li(7,ii);x_li(8,ii);x_li(9,ii)];
    ea_d = [x_li(1,ii);x_li(2,ii);x_li(3,ii)];
    vR_d = reshape(expm(hatMap(ea_d)),[9,1]);
    pfd = reshape(Rground * pf34,[12,1]);
    Xd_(:,ii) = [pc_d;dpc_d;vR_d;wb_d;pfd];
    Ud_(:,ii) = f_li(:,ii);
end
end

%% Tool functions
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
function R=rotsb(theta)% constructing rotation matrix
% R=rotx(theta(1))*roty(theta(2))*rotz(theta(3));
R=rotz(theta(3))*roty(theta(2))*rotx(theta(1));

end

function s=Skew(in)
s = [0 -in(3) in(2);
    in(3) 0 -in(1);
    -in(2) in(1) 0];
end