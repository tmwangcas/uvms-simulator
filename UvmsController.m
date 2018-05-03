classdef UvmsController
    
    properties
        uvms_dynamics;
        vf, vf_ee;
        err2_init, err_ee2_init;
    end
    
    properties(Constant)
        % controller parameters        
        alpha1 = 0.5;
        alpha2 = 0.5;
        beta   = 1;
        ks     = 500;
        
%         alpha1 = 5;
%         alpha2 = 5;
%         beta   = 1;
%         ks     = 500;
    end
    
    methods
        function obj = UvmsController(uvms_dynamics)
            if nargin == 1
                obj.uvms_dynamics = uvms_dynamics;
                obj.vf    = zeros(9,1);
                obj.vf_ee = zeros(3,1);
                obj.err2_init    = zeros(9,1);
                obj.err_ee2_init = zeros(3,1);
            else
                error('Not enough parameters');
            end
        end
        
        %% Calculate Control Inputs
        function obj = GetControlInputs(obj)
            J_veh_vel = obj.uvms_dynamics.uvms_kinematics.get_J_veh_vel(); % Jacobian matrix for vehicle
            err2  = obj.GetTrackingError();
            obj = obj.get_vf();
            
            if obj.uvms_dynamics.uvms_kinematics.t == 0
                obj.err2_init = err2;
            end
            
            tau     = (obj.ks+1)*err2 - (obj.ks+1)*obj.err2_init + obj.vf; % 9 DOF control forces in inertial frame
            tau_v_I = tau(1:6);
            tau_m   = tau(7:9);
            tau_v_0 = J_veh_vel'*tau_v_I;
            obj.uvms_dynamics.tau_c = [tau_v_0; tau_m]; % 9 DOF generalized control forces
            
            % TEST
%             obj.uvms_dynamics.tau_c = zeros(9,1);
%             obj.uvms_dynamics.tau_c = [0; 0; 0; 0; 0; 0; 0.5; 0; 0];
        end
        
        function obj = GetControlInputs_ee(obj)
            J_dk_lin = obj.uvms_dynamics.uvms_kinematics.get_J_dk_lin();
            err_ee2  = obj.GetEndEffectorTrackingError();
            obj = obj.get_vf_ee();
            
            if obj.uvms_dynamics.uvms_kinematics.t == 0
                obj.err_ee2_init = err_ee2;
            end
            
            tau_ee = (obj.ks+1)*err_ee2 - (obj.ks+1)*obj.err_ee2_init + obj.vf_ee; % 3 DOF control forces in end-effector frame
            obj.uvms_dynamics.tau_c = J_dk_lin'*tau_ee; % 9 DOF generalized control forces
            
%             obj.uvms_dynamics.tau_c = obj.uvms_dynamics.tau_c + normrnd(0,0.01,[9,1]); % Gaussian noise
            
            % TEST
%             obj.uvms_dynamics.tau_c = zeros(9,1);
%             obj.uvms_dynamics.tau_c = [0; 0; 0; 0; 0; 0; 0.5; 0; 0];
        end
        
        % calculate vf
        function obj = get_vf(obj)
            err2 = obj.GetTrackingError();
            dvf = (obj.ks+1)*obj.alpha2*err2 + obj.beta*sign(err2);
            
            if obj.uvms_dynamics.uvms_kinematics.t == 0
                obj.vf = zeros(9,1);
            else
                obj.vf = obj.vf + dvf*obj.uvms_dynamics.uvms_kinematics.dt;
            end
        end
        
        % calculate vf_ee
        function obj = get_vf_ee(obj)
            err_ee2 = obj.GetEndEffectorTrackingError();
            dvf_ee = (obj.ks+1)*obj.alpha2*err_ee2 + obj.beta*sign(err_ee2);
            
            if obj.uvms_dynamics.uvms_kinematics.t == 0
                obj.vf_ee = zeros(3,1);
            else
                obj.vf_ee = obj.vf_ee + dvf_ee*obj.uvms_dynamics.uvms_kinematics.dt;
            end
        end
        
        %% Calculate Tracking Error
        % tracking error of generalized coordinates and their derivative
        function err2 = GetTrackingError(obj)
            [eta_d, deta_d] = obj.uvms_dynamics.uvms_kinematics.GetDesiredTrajectory();
            err1  = eta_d - obj.uvms_dynamics.uvms_kinematics.eta;
            derr1 = deta_d - obj.uvms_dynamics.uvms_kinematics.deta;
            err2  = derr1 + obj.alpha1*err1;
        end
        
        % tracking error of end-effector position and linear velocity
        function err_ee2 = GetEndEffectorTrackingError(obj)
            [p_ee_d, dp_ee_d] = obj.uvms_dynamics.uvms_kinematics.GetEndEffectorDesiredTrajectory();
            [p_ee, dp_ee]     = obj.uvms_dynamics.uvms_kinematics.DirectKinematics();
            err_ee1  = p_ee_d - p_ee;
            derr_ee1 = dp_ee_d - dp_ee;
            err_ee2  = derr_ee1 + obj.alpha1*err_ee1;
        end
    end
end