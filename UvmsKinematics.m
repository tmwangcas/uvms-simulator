classdef UvmsKinematics
    
    properties
        p, theta, q;
        dp, dtheta;
        v, omega, dq;
        dv, domega, ddq;
        eta, deta, zeta, dzeta;
        eta_init, deta_init;
        p_ee_init, dp_ee_init;
        t, dt;
    end
    
    properties(Constant)        
        %% Robot Geometry
        % dimension of vehicle
        Lx = 0.9;   % length
        Ly = 0.79;  % width
        Lz = 0.249; % height
        Ly_notch = 0.182; % width of notch on vehicle

        % length of link
        L1 = 0.17;
        L2 = 0.29;
        L3 = 0.29;
        L1_horizontal = 0.1573; % length of horizontal cylinder of link 1
        
        % radius of link
        r1 = 0.042;
        r2 = 19.8*1e-3;
        r3 = 7.9*1e-3;
        
        %% Coordinate System
        % vector of coordinate axis
        z_1_0 = [0; 0; 1];
        y_2_1 = [0; 1; 0];
        
        z_1_1 = [0; 0; 1];
        z_2_2 = [0; 0; 1];
        z_3_3 = [0; 0; 1];
        
        % distance between coordinate origins
        do_1_0_x = 0.044;
        do_1_0_z = 0.01;
        do_3_2_x = 0.07;
        
        % radius of end-effector desired trajectory
        r_ee_traj = 0.24;
    end
    
    methods
        function obj = UvmsKinematics(dt)
            if nargin == 1
                % initial state
                obj.p      = zeros(3,1);
                obj.theta  = zeros(3,1);
                obj.q      = [0; pi/2; 0.32];
                
                obj.dp     = zeros(3,1);
                obj.dtheta = zeros(3,1);
                
                obj.v      = zeros(3,1);
                obj.omega  = zeros(3,1);
                obj.dq     = zeros(3,1);
                
                obj.dv     = zeros(3,1);
                obj.domega = zeros(3,1);
                obj.ddq    = zeros(3,1);
                
                obj.eta    = [ obj.p;  obj.theta;   obj.q];
                obj.deta   = [obj.dp; obj.dtheta;  obj.dq];
                obj.zeta   = [ obj.v;  obj.omega;  obj.dq];
                obj.dzeta  = [obj.dv; obj.domega; obj.ddq];
                
                obj.eta_init  = obj.eta;
                obj.deta_init = obj.deta;
                
                [obj.p_ee_init, obj.dp_ee_init] = obj.DirectKinematics();
                
                obj.dt     = dt;
                obj.t      = 0;
            else
                error('Not enough parameters');
            end
        end
        
        %% Update Robot State
        function obj = UpdateKinematics(obj)
            % transformation matrix
            R_0_I         = obj.get_R_0_I();
            J_veh_ang_vel = obj.get_J_veh_ang_vel();
            
            % update kinematics
            obj.dv     = obj.dzeta(1:3);
            obj.domega = obj.dzeta(4:6);
            obj.ddq    = obj.dzeta(7:9);
            
            obj.dp     = R_0_I*obj.v;
            obj.dtheta = J_veh_ang_vel*obj.omega;
            
            obj.p     = obj.p     + obj.dp*obj.dt;
            obj.theta = obj.theta + obj.dtheta*obj.dt;
            obj.q     = obj.q     + obj.dq*obj.dt;
            
            % range of q3 is [0.32, 0.58]
%             if obj.q(3) < 0.32
%                 obj.q(3) = 0.32;
%                 obj.dq(3) = 0;
%             elseif obj.q(3) > 0.58
%                 obj.q(3) = 0.58;
%                 obj.dq(3) = 0;
%             end
                
            obj.v     = obj.v     + obj.dv*obj.dt;          
            obj.omega = obj.omega + obj.domega*obj.dt;
            obj.dq    = obj.dq    + obj.ddq*obj.dt;
            
            obj.eta   = [ obj.p;  obj.theta;  obj.q];
            obj.deta  = [obj.dp; obj.dtheta; obj.dq];
            obj.zeta  = [ obj.v;  obj.omega; obj.dq];
            
            % update time
            obj.t = obj.t + obj.dt;
        end
        
        %% Get Desired Trajectory
        % desired value of generalized ccordinates and their derivative
        function [eta_d, deta_d] = GetDesiredTrajectory(obj)
            % case 1
            eta_d  = [0; 0; 0; 0; 0; 0; 0; pi/4; 0.58];
            deta_d = zeros(9,1);
            
            % case 2
            tt = 0.5;
            if obj.t >= 0 && obj.t < tt
                eta_d  = [0; 0; 0; 0; 0; 0; 0; pi/4; 0.45];
                deta_d = [0; 0; 0; 0; 0; 0; 0; -pi/4/tt; 0.13/tt];
            elseif obj.t >= tt && obj.t < tt*2
                eta_d  = [0; 0; 0; 0; 0; 0; pi/4; pi/4; 0.58];
                deta_d = [0; 0; 0; 0; 0; 0; pi/4/tt; 0; 0.13/tt];
            elseif obj.t >= tt*2 && obj.t < tt*3
                eta_d  = [0; 0; 0; 0; 0; 0; pi/4; pi/2; 0.45];
                deta_d = [0; 0; 0; 0; 0; 0; 0; pi/4/tt; -0.13/tt];
            else
                eta_d  = [0; 0; 0; 0; 0; 0; 0; pi/2; 0.32];
                deta_d = [0; 0; 0; 0; 0; 0; -pi/4/tt; 0; -0.13/tt];
            end
            
            % TEST
%             eta_d  = obj.eta_init;
%             deta_d = obj.deta_init;
        end
        
        % desired value of position and linear velocity of end-effector
        function [p_ee_d, dp_ee_d] = GetEndEffectorDesiredTrajectory(obj)            
            omega_ee_traj = 2*pi/20; % angular velocity, period T = 10s
            p_ee_d  = obj.p_ee_init + [0; 0; obj.r_ee_traj] + [0; obj.r_ee_traj*cos(omega_ee_traj*obj.t-pi/2); obj.r_ee_traj*sin(omega_ee_traj*obj.t-pi/2)];
            dp_ee_d = obj.dp_ee_init + [0; -omega_ee_traj*obj.r_ee_traj*sin(omega_ee_traj*obj.t-pi/2); omega_ee_traj*obj.r_ee_traj*cos(omega_ee_traj*obj.t-pi/2)];
        end
        
        %% Direct Kinematics        
        % Direct Kinematics (DK)
        % calculate end-effector position expressed in inertial frame from generalized coordinates
        % p_ee = po_3_I
        % Direct Differential Kinematics (DDK)
        % calculate end-effector linear velocity expressed in inertial frame from generalized velocity
        % dp_ee = dpo_3_I
        function [p_ee, dp_ee] = DirectKinematics(obj)
            R_0_I = obj.get_R_0_I();
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            [po_1_0, po_2_1, po_3_2] = obj.GetCoordinateOrigin();
            J_dk_lin = obj.get_J_dk_lin();
            
            p_ee = R_0_I*(R_1_0*(R_2_1*po_3_2 + po_2_1) + po_1_0) + obj.p;
            dp_ee = J_dk_lin*obj.zeta;
        end
        
        %% Inverse Kinematics
%         % Inverse Kinematics (IK)
%         % calculate generalized coordinates from end-effector position expressed in inertial frame
%         % p_ee_d = po_3_I_d
%         % Inverse Differential Kinematics (IDK)
%         % calculate generalized velocity from end-effector linear velocity expressed in inertial frame
%         % dp_ee_d = dpo_3_I_d
%         function [eta_d, zeta_d] = InverseKinematics(obj, p_ee_d, dp_ee_d)
%             J_k_pos = obj.get_J_k_pos();
%             J_dk_lin = obj.get_J_dk_lin();
%             
%             eta_d = pinv(J_k_pos)*p_ee_d;
%             zeta_d = pinv(J_dk_lin)*dp_ee_d;
%         end
        
        %% UVMS Geometry  
        % Vehicle Geometry
        % position vector of vehicle vertex expressed in inertial frame
        function vehicle_vertex_I = GetVehicleGeometry(obj)
            R_0_I = obj.get_R_0_I();
            
            vehicle_vertex_orientation = [1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1;
                                          1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1;
                                          1  1  1  1 -1 -1 -1 -1  1  1  1  1  1  1  1  1];
            vehicle_vertex_distance = [[obj.Lx/2; obj.Ly/2; obj.Lz/2]*ones(1,8) ...
                [obj.Lx/2; obj.Ly_notch/2; obj.Lz/2]*ones(1,4) ...
                [obj.Lx/2; obj.Ly_notch/2; obj.do_1_0_z]*ones(1,4)];
            vehicle_vertex_0 = vehicle_vertex_orientation.*vehicle_vertex_distance;
            
            vehicle_vertex_I = R_0_I*vehicle_vertex_0 + obj.p*ones(1,16);
        end
        
        % Manipulator Geometry
        % position vector of several specific points on manipulator expressed in inertial frame
        function [po_1_I, po_2_I, po_3_I, p_L1_left_I, p_L1_right_I, p_L2_start_I, p_L2_end_I] = GetManipulatorGeometry(obj)
            R_0_I = obj.get_R_0_I();
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            [po_1_0, po_2_1, po_3_2] = obj.GetCoordinateOrigin();
            
            po_1_I = R_0_I*po_1_0 + obj.p;
            po_2_I = R_0_I*(R_1_0*po_2_1 + po_1_0) + obj.p;
            po_3_I = R_0_I*(R_1_0*(R_2_1*po_3_2 + po_2_1) + po_1_0) + obj.p;
            
            p_L1_left_2  = [0;  obj.L1_horizontal/2; 0];
            p_L1_right_2 = [0; -obj.L1_horizontal/2; 0];
            p_L1_left_I  = R_0_I*(R_1_0*(R_2_1*p_L1_left_2 + po_2_1) + po_1_0) + obj.p;
            p_L1_right_I = R_0_I*(R_1_0*(R_2_1*p_L1_right_2 + po_2_1) + po_1_0) + obj.p;

            p_L2_start_2 = [-obj.do_3_2_x; 0; 0];
            p_L2_end_2   = [-obj.do_3_2_x; 0; obj.L2];
            p_L2_start_I = R_0_I*(R_1_0*(R_2_1*p_L2_start_2 + po_2_1) + po_1_0) + obj.p;
            p_L2_end_I   = R_0_I*(R_1_0*(R_2_1*p_L2_end_2 + po_2_1) + po_1_0) + obj.p;
        end
        
        %% Center of Mass
        % position vector of C.M. of link 0 1 2 3
        function [pc_0_0, pc_1_1, pc_2_2, pc_3_3] = GetCenterOfMass(obj)
            pc_0_0 = [0; 0; 0];
            pc_1_1 = [0; 0; 0.12];
            pc_2_2 = [-obj.do_3_2_x; 0; obj.L2/2];
            pc_3_3 = [0; 0; -obj.L3/2];
        end
        
        %% Distance between C.M.
        % position vector of C.M. of link 1 2 3 with respect to C.M. of link 0 expressed in frame 0
        function [dc_1_0, dc_2_0, dc_3_0] = GetDistanceBetweenCM(obj)
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            [po_1_0, po_2_1, po_3_2] = obj.GetCoordinateOrigin();
            [pc_0_0, pc_1_1, pc_2_2, pc_3_3] = obj.GetCenterOfMass();
            
            dc_1_0 = R_1_0*pc_1_1 + po_1_0 - pc_0_0;
            dc_2_0 = R_1_0*(R_2_1*pc_2_2 + po_2_1) + po_1_0 - pc_0_0;
            dc_3_0 = R_1_0*(R_2_1*(pc_3_3 + po_3_2) + po_2_1) + po_1_0 - pc_0_0;
        end
        
        %% Coordinate Origin
        % position vector of origin of frame 1 2 3 expressed in frame 0 1 2
        function [po_1_0, po_2_1, po_3_2] = GetCoordinateOrigin(obj)
            po_1_0 = [-obj.do_1_0_x; 0; obj.do_1_0_z];
            po_2_1 = [0; 0; obj.L2];
            po_3_2 = [-obj.do_3_2_x; 0; obj.q(3)];
        end
        
        %% Rotation Matrix
        % rotation matrix from frame 0 to inertial frame (axis rotation sequence is XYZ)
        function R_0_I = get_R_0_I(obj) 
            s1 = sin(obj.theta(1));  c1 = cos(obj.theta(1));
            s2 = sin(obj.theta(2));  c2 = cos(obj.theta(2));
            s3 = sin(obj.theta(3));  c3 = cos(obj.theta(3));
                  
            R_0_I = [c2*c3  -c1*s3+s1*s2*c3   s1*s3+c1*s2*c3                  
                     c2*s3   c1*c3+s1*s2*s3  -s1*c3+c1*s2*s3
                       -s2            s1*c2            c1*c2];
        end
        
        % rotation matrix from frame 1 to frame 0
        function R_1_0 = get_R_1_0(obj)
            R_1_0 = [cos(obj.q(1)) -sin(obj.q(1))  0
                     sin(obj.q(1))  cos(obj.q(1))  0
                                 0              0  1];
        end
        
        % rotation matrix from frame 2 to frame 1
        function R_2_1 = get_R_2_1(obj)
            R_2_1 = [ cos(obj.q(2))  0  sin(obj.q(2))
                                  0  1              0
                     -sin(obj.q(2))  0  cos(obj.q(2))];
        end
        
        % rotation matrix from frame 2 to frame 0
        function R_2_0 = get_R_2_0(obj)
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            
            R_2_0 = R_1_0*R_2_1;
        end
        
        % rotation matrix from frame 3 to frame 0
        function R_3_0 = get_R_3_0(obj)
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            R_3_2 = eye(3);
            
            R_3_0 = R_1_0*R_2_1*R_3_2;
        end
        
        %% Jacobian Matrix        
        % Jacobian matrix to transform body-fixed angular velocity to derivative of Euler angles
        function J_veh_ang_vel = get_J_veh_ang_vel(obj)
            s1 = sin(obj.theta(1));  c1 = cos(obj.theta(1));
            c2 = cos(obj.theta(2));  t2 = tan(obj.theta(2));
            
            J_veh_ang_vel = [1  s1*t2  c1*t2
                             0     c1    -s1
                             0  s1/c2  c1/c2];
            
            % J_veh_ang_vel is singular for obj.theta(2) == (2k+1)*pi/2
            if obj.theta(2) == pi/2 || obj.theta(2) == -pi/2
                J_veh_ang_vel = eye(3);
            end
        end
        
        % Jacobian matrix to transform body-fixed linear & angular velocities to derivative of position & Euler angles
        function J_veh_vel = get_J_veh_vel(obj)
            R_0_I         = obj.get_R_0_I();
            J_veh_ang_vel = obj.get_J_veh_ang_vel();
            
            J_veh_vel = [     R_0_I     zeros(3,3)
                         zeros(3,3)  J_veh_ang_vel];
        end
        
        % Direct Differential Kinematics (DDK): J_dk_lin
        % Jacobian matrix to transform generalized velocity to end-effector linear velocity expressed in inertial frame
        % Inverse Differential Kinematics (IDK): pinv(J_dk_lin)
        % Jacobian matrix to transform end-effector linear velocity expressed in inertial frame to generalized velocity
        function J_dk_lin = get_J_dk_lin(obj)
            R_0_I = obj.get_R_0_I();
            R_1_0 = obj.get_R_1_0();
            R_2_1 = obj.get_R_2_1();
            [po_1_0, po_2_1, po_3_2] = obj.GetCoordinateOrigin();
            
            po_3_0 = R_1_0*(R_2_1*po_3_2 + po_2_1) + po_1_0;
            J_dpo_3_0 = [obj.S(obj.z_1_0)*R_1_0*(R_2_1*po_3_2+po_2_1) R_1_0*obj.S(obj.y_2_1)*R_2_1*po_3_2 R_1_0*R_2_1*obj.z_2_2];
            J_dk_lin = [R_0_I -obj.S(R_0_I*po_3_0)*R_0_I R_0_I*J_dpo_3_0];
        end
    end
    
    methods(Static)
        % matrix operator performing cross between two vectors
        function out = S(x)
            out = [    0 -x(3)  x(2)
                    x(3)     0 -x(1)
                   -x(2)  x(1)     0];
        end
    end  
end