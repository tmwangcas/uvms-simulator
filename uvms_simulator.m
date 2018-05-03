clc
clearvars
close all

%% Simulation Parameters
t_f   = 5;        % final simulation time [s]
dt    = 0.001;     % sampling time [s]
t_vec = 0:dt:t_f;  % time vector [s]
t_num = length(t_vec); % number of samples

%% Definition of Variables
p     = zeros(3,t_num);
theta = zeros(3,t_num);
q     = zeros(3,t_num);

dp     = zeros(3,t_num);
dtheta = zeros(3,t_num);

v     = zeros(3,t_num);
omega = zeros(3,t_num);
dq    = zeros(3,t_num);

dv     = zeros(3,t_num);
domega = zeros(3,t_num);
ddq    = zeros(3,t_num);

eta   = zeros(9,t_num); % generalized coordinates
deta  = zeros(9,t_num);
zeta  = zeros(9,t_num); % generalized velocity
dzeta = zeros(9,t_num); % generalized acceleration

p_ee  = zeros(3,t_num);
dp_ee = zeros(3,t_num);

p_ee_d  = zeros(3,t_num);
dp_ee_d = zeros(3,t_num);

po_1_I = zeros(3,t_num);
po_2_I = zeros(3,t_num);
po_3_I = zeros(3,t_num);

p_L2_left_I  = zeros(3,t_num);
p_L2_right_I = zeros(3,t_num);
p_L2_start_I = zeros(3,t_num);
p_L2_end_I   = zeros(3,t_num);

vf_I = zeros(3,t_num);
af_I = zeros(3,t_num);

tau_th = zeros(8,t_num); % thruster control torques
tau_v  = zeros(6,t_num); % vehicle control forces in vehicle frame
tau_m  = zeros(3,t_num); % manipulator control forces
tau_c  = zeros(9,t_num);

%% Kinematics
uvms_kinematics = UvmsKinematics(dt);

%% Dynamics
uvms_dynamics = UvmsDynamics(uvms_kinematics);

%% Controller
uvms_controller = UvmsController(uvms_dynamics);

%% Graphics
uvms_graphics = UvmsGraphics(uvms_dynamics);

%% Display Robot Movement
f1 = figure;
set(gcf,'outerposition',get(0,'screensize'))

%% Main Loop
for i = 1:t_num
    %% Current Time
    t_cur = (i-1)*dt;
    
    %% Fluid Velocity & Acceleration
    vf_I(:,i) = zeros(3,1);
    af_I(:,i) = zeros(3,1);
    
%     A = 0.3; % amplitude
%     w = 2*pi/5; % angular velocity
%     vf_I(:,i) = [0; A*sin(w*t_cur); 0];
%     af_I(:,i) = [0; A*w*cos(w*t_cur); 0];
    
    %% Controller
%     uvms_controller = uvms_controller.GetControlInputs();
    uvms_controller = uvms_controller.GetControlInputs_ee();
    
    %% Direct Dynamics
    uvms_controller.uvms_dynamics = uvms_controller.uvms_dynamics.DirectDynamics(vf_I(:,i), af_I(:,i));
    
    %% Get Current State Variable
    p(:,i)     = uvms_controller.uvms_dynamics.uvms_kinematics.p;
    theta(:,i) = uvms_controller.uvms_dynamics.uvms_kinematics.theta;
    q(:,i)     = uvms_controller.uvms_dynamics.uvms_kinematics.q;
    
    v(:,i)     = uvms_controller.uvms_dynamics.uvms_kinematics.v;
    omega(:,i) = uvms_controller.uvms_dynamics.uvms_kinematics.omega;
    dq(:,i)    = uvms_controller.uvms_dynamics.uvms_kinematics.dq;
    
    dv(:,i)     = uvms_controller.uvms_dynamics.uvms_kinematics.dv;
    domega(:,i) = uvms_controller.uvms_dynamics.uvms_kinematics.domega;
    ddq(:,i)    = uvms_controller.uvms_dynamics.uvms_kinematics.ddq;
    
    eta(:,i)   = uvms_controller.uvms_dynamics.uvms_kinematics.eta;
    zeta(:,i)  = uvms_controller.uvms_dynamics.uvms_kinematics.zeta;
    dzeta(:,i) = uvms_controller.uvms_dynamics.uvms_kinematics.dzeta;
    
    tau_c(:,i) = uvms_controller.uvms_dynamics.tau_c;
    
    [p_ee(:,i), dp_ee(:,i)]     = uvms_controller.uvms_dynamics.uvms_kinematics.DirectKinematics();
    [p_ee_d(:,i), dp_ee_d(:,i)] = uvms_controller.uvms_dynamics.uvms_kinematics.GetEndEffectorDesiredTrajectory();
    
    p_ee_init  = uvms_controller.uvms_dynamics.uvms_kinematics.p_ee_init;
    dp_ee_init = uvms_controller.uvms_dynamics.uvms_kinematics.dp_ee_init;
    
    %% Update Graphics Variables
    uvms_graphics = uvms_graphics.UpdateUvmsDynamics(uvms_controller.uvms_dynamics);
    
    %% Display Robot Movement
    if mod(t_cur,0.01) == 0
        figure(f1);
        
        uvms_graphics.DrawEndEffectorDesiredTrajectory();
        hold on
        uvms_graphics.DrawVehicle();
        uvms_graphics.DrawManipulator();
        hold off
        
        xlabel('X position/m')
        ylabel('Y position/m')
        zlabel('Z position/m')
%         axis([-0.8 0.8 -0.8 0.8 -0.3 0.9])
        axis([-1.6 1.6 -1.6 1.6 -0.6 1.8])
%         axis equal
        grid on
        title('3D Trajectory')
        view(60,10)
        drawnow
    end
    
    %% Update Kinematics
    uvms_controller.uvms_dynamics.uvms_kinematics = uvms_controller.uvms_dynamics.uvms_kinematics.UpdateKinematics();
    
    %% Count Loop
    disp(i)
end

%% Plot 
% draw end-effector trajectory
uvms_graphics.DrawEndEffectorTracjectory(p_ee, p_ee_d);

% plot end-effector trajectory
uvms_graphics.PlotEndEffectorTracjectory(t_vec, p_ee, p_ee_d);

% plot generalized coordinate
uvms_graphics.PlotGeneralizedCoordinate(t_vec, eta);

% plot generalized velocity
uvms_graphics.PlotGeneralizedVelocity(t_vec, zeta);

% plot control forces
uvms_graphics.PlotControlForces(t_vec, tau_c);

