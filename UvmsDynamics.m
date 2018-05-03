classdef UvmsDynamics
    
    properties
        uvms_kinematics;
        tau_c;
    end
    
    properties(Constant)        
        %% Robot System Parameters
        % mass of link
        m0 = 60; 
        m1 = 4.2;
        m2 = 1;
        m3 = 1;

        % inertia matrix of link
        I_0_0 = diag([3.368 3.553 6.244]);
        I_1_1 = diag([19.11*1e-3 15.56*1e-3  7.54*1e-3]);
        I_2_2 = diag([13.97*1e-3 24.66*1e-3 12.19*1e-3]);
        I_3_3 = diag([13.97*1e-3 24.66*1e-3 12.19*1e-3]);
        
        %% Gravitational Acceleration
        g_I = [0; 0; -9.81];       
        
        %% Hydrodynamic Effects
        % fluid density
        rho = 1000;
        
        % drag coefficients
        Cd0 = 0.47; % sphere
        Cd1 = 1.1;  % cylinder
        Cd2 = 1.1;  % cylinder
        Cd3 = 1.1;  % cylinder
    end
    
    methods
        function obj = UvmsDynamics(uvms_kinematics)
            if nargin == 1
                obj.uvms_kinematics = uvms_kinematics;
                obj.tau_c = zeros(9,1);
            else
                error('Not enough parameters');
            end
        end
        
        %% Direct Dynamics
        function obj = DirectDynamics(obj, vf_I, af_I)   
            [M, C, Fe] = obj.GetModelParameters(vf_I, af_I);
            obj.uvms_kinematics.dzeta = M\(obj.tau_c - C - Fe);
        end
        
        %% Model Parameters
        function [M, C, Fe] = GetModelParameters(obj, vf_I, af_I)
            % generalized velocity
            v     = obj.uvms_kinematics.v;
            omega = obj.uvms_kinematics.omega;
            dq1   = obj.uvms_kinematics.dq(1);
            dq2   = obj.uvms_kinematics.dq(2);
            dq3   = obj.uvms_kinematics.dq(3);
            zeta  = obj.uvms_kinematics.zeta;
            
            % rotation matrix
            R_0_I = obj.uvms_kinematics.get_R_0_I();
            R_1_0 = obj.uvms_kinematics.get_R_1_0();
            R_2_1 = obj.uvms_kinematics.get_R_2_1();
            R_2_0 = obj.uvms_kinematics.get_R_2_0();
            R_3_0 = obj.uvms_kinematics.get_R_3_0();
            
            % coordinate origin
            [po_1_0, po_2_1, po_3_2] = obj.uvms_kinematics.GetCoordinateOrigin();
            
            % position vector of center of mass            
            [pc_0_0, pc_1_1, pc_2_2, pc_3_3] = obj.uvms_kinematics.GetCenterOfMass();
            
            % distance vector between C.M.
            [dc_1_0, dc_2_0, dc_3_0] = obj.uvms_kinematics.GetDistanceBetweenCM();
            
            % radius of vehicle
            r0 = (obj.uvms_kinematics.Lx*obj.uvms_kinematics.Ly*obj.uvms_kinematics.Lz/8)^(1/3);
            
            % volume of link
            [V0, V1, V2, V3] = obj.GetVolume();
            
            % inertia matrix of link                        
            [I_1_0, I_2_0, I_3_0] = obj.GetInertiaMatrix();

            % added mass matrix
            [Ia_0_0, Ia_1_0, Ia_2_0, Ia_3_0] = obj.GetAddedMassMatrix();
            
            % gravitational acceleration
            g_0 = R_0_I'*obj.g_I;
            
            % fluid velocity & acceleration
            vf_0 = R_0_I'*vf_I;
            af_0 = R_0_I'*af_I;
            
            %% Coefficients of Derivative of Position Vector 
            % dp
            A_dq1_1 = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*pc_1_1;

            A_dq1_2 = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*(R_2_1*pc_2_2 + po_2_1);
            A_dq2_2 = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*pc_2_2;

            A_dq1_3 = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*(R_2_1*(pc_3_3 + po_3_2) + po_2_1);
            A_dq2_3 = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*(pc_3_3 + po_3_2);
            A_dq3_3 = R_1_0*R_2_1*obj.uvms_kinematics.z_2_2;

            % ddp
            B_ddq1_1   = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*pc_1_1;
            B_dq1sqr_1 = obj.S(obj.uvms_kinematics.z_1_0)^2*R_1_0*pc_1_1;

            B_ddq1_2   = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*(R_2_1*pc_2_2 + po_2_1);
            B_ddq2_2   = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*pc_2_2;
            B_dq1sqr_2 = obj.S(obj.uvms_kinematics.z_1_0)^2*R_1_0*(R_2_1*pc_2_2 + po_2_1);
            B_dq2sqr_2 = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)^2*R_2_1*pc_2_2;
            B_dq1dq2_2 = 2*obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*pc_2_2;

            B_ddq1_3   = obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*(R_2_1*(pc_3_3 + po_3_2) + po_2_1);
            B_ddq2_3   = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*(pc_3_3 + po_3_2);
            B_ddq3_3   = R_1_0*R_2_1*obj.uvms_kinematics.z_2_2;
            B_dq1sqr_3 = obj.S(obj.uvms_kinematics.z_1_0)^2*R_1_0*(R_2_1*(pc_3_3 + po_3_2) + po_2_1);
            B_dq2sqr_3 = R_1_0*obj.S(obj.uvms_kinematics.y_2_1)^2*R_2_1*(pc_3_3 + po_3_2);
            B_dq1dq2_3 = 2*obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*(pc_3_3 + po_3_2);
            B_dq1dq3_3 = 2*obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*R_2_1*obj.uvms_kinematics.z_2_2;
            B_dq2dq3_3 = 2*R_1_0*obj.S(obj.uvms_kinematics.y_2_1)*R_2_1*obj.uvms_kinematics.z_2_2;
            
            %% Jacobian Matrix
            % angular velocity
            J_omega_0 = [zeros(3,3) eye(3) zeros(3,3)];
            J_omega_1 = [zeros(3,3) eye(3) obj.uvms_kinematics.z_1_0 zeros(3,2)];
            J_omega_2 = [zeros(3,3) eye(3) obj.uvms_kinematics.z_1_0 R_1_0*obj.uvms_kinematics.y_2_1 zeros(3,1)];
            J_omega_3 = J_omega_2;

            % linear velocity
            J_v_0 = [eye(3) zeros(3,3) zeros(3,3)];
            J_v_1 = [eye(3) -obj.S(dc_1_0) A_dq1_1 zeros(3,2)];
            J_v_2 = [eye(3) -obj.S(dc_2_0) A_dq1_2 A_dq2_2 zeros(3,1)];
            J_v_3 = [eye(3) -obj.S(dc_3_0) A_dq1_3 A_dq2_3 A_dq3_3];

            % angular acceleration
            J_alpha_dzeta_0 = [zeros(3,3) eye(3) zeros(3,3)];
            J_alpha_dzeta_1 = [zeros(3,3) eye(3) obj.uvms_kinematics.z_1_0 zeros(3,2)];
            J_alpha_dzeta_2 = [zeros(3,3) eye(3) obj.uvms_kinematics.z_1_0 R_1_0*obj.uvms_kinematics.y_2_1 zeros(3,1)];
            J_alpha_dzeta_3 = J_alpha_dzeta_2;

            J_alpha_zeta_0 = zeros(3,9);
            J_alpha_zeta_1 = [zeros(3,3) -0.5*obj.S(obj.uvms_kinematics.z_1_0)*dq1 -0.5*obj.S(obj.uvms_kinematics.z_1_0)*omega zeros(3,2)];
            J_alpha_zeta_2 = [zeros(3,3) -0.5*(obj.S(obj.uvms_kinematics.z_1_0)*dq1+obj.S(R_1_0*obj.uvms_kinematics.y_2_1)*dq2) ...
                -0.5*obj.S(obj.uvms_kinematics.z_1_0)*(omega-R_1_0*obj.uvms_kinematics.y_2_1*dq2) ...
                -0.5*(obj.S(R_1_0*obj.uvms_kinematics.y_2_1)*omega-obj.S(obj.uvms_kinematics.z_1_0)*R_1_0*obj.uvms_kinematics.y_2_1*dq1) zeros(3,1)];
            J_alpha_zeta_3 = J_alpha_zeta_2;

            % linear acceleration
            J_a_dzeta_0 = [eye(3) zeros(3,3) zeros(3,3)];
            J_a_dzeta_1 = [eye(3) -obj.S(dc_1_0) B_ddq1_1 zeros(3,2)];
            J_a_dzeta_2 = [eye(3) -obj.S(dc_2_0) B_ddq1_2 B_ddq2_2 zeros(3,1)];
            J_a_dzeta_3 = [eye(3) -obj.S(dc_3_0) B_ddq1_3 B_ddq2_3 B_ddq3_3];

            J_a_zeta_0 = [0.5*obj.S(omega) -0.5*obj.S(v) zeros(3,3)];
            J_a_zeta_1 = [0.5*obj.S(omega) -(0.5*obj.S(v)+obj.S(A_dq1_1)*dq1) ...
                -(obj.S(A_dq1_1)*omega-B_dq1sqr_1*dq1) zeros(3,2)];
            J_a_zeta_2 = [0.5*obj.S(omega) -(0.5*obj.S(v)+obj.S(A_dq1_2)*dq1+obj.S(A_dq2_2)*dq2) ...
                -(obj.S(A_dq1_2)*omega-B_dq1sqr_2*dq1-0.5*B_dq1dq2_2*dq2) ...
                -(obj.S(A_dq2_2)*omega-B_dq2sqr_2*dq2-0.5*B_dq1dq2_2*dq1) zeros(3,1)];
            J_a_zeta_3 = [0.5*obj.S(omega) -(0.5*obj.S(v)+obj.S(A_dq1_3)*dq1+obj.S(A_dq2_3)*dq2+obj.S(A_dq3_3)*dq3) ...
                -(obj.S(A_dq1_3)*omega-B_dq1sqr_3*dq1-0.5*B_dq1dq2_3*dq2-0.5*B_dq1dq3_3*dq3) ...
                -(obj.S(A_dq2_3)*omega-B_dq2sqr_3*dq2-0.5*B_dq1dq2_3*dq1-0.5*B_dq2dq3_3*dq3) ...
                -(obj.S(A_dq3_3)*omega-0.5*B_dq1dq3_3*dq1-0.5*B_dq2dq3_3*dq2)];
            
            %% Model Parameters
            % inertia matrix (including added mass)
            M0 = [J_v_0; J_omega_0]'*([obj.m0*J_a_dzeta_0; obj.I_0_0*J_alpha_dzeta_0] + Ia_0_0*[J_a_dzeta_0; J_alpha_dzeta_0]);
            M1 = [J_v_1; J_omega_1]'*([obj.m1*J_a_dzeta_1; I_1_0*J_alpha_dzeta_1] + Ia_1_0*[J_a_dzeta_1; J_alpha_dzeta_1]);
            M2 = [J_v_2; J_omega_2]'*([obj.m2*J_a_dzeta_2; I_2_0*J_alpha_dzeta_2] + Ia_2_0*[J_a_dzeta_2; J_alpha_dzeta_2]);
            M3 = [J_v_3; J_omega_3]'*([obj.m3*J_a_dzeta_3; I_3_0*J_alpha_dzeta_3] + Ia_3_0*[J_a_dzeta_3; J_alpha_dzeta_3]);

            M = M0 + M1 + M2 + M3;

            % vector of Coriolis and centripetal terms (including added mass)
            C0 = [J_v_0; J_omega_0]'*([obj.m0*J_a_zeta_0; obj.I_0_0*J_alpha_zeta_0+obj.S(J_omega_0*zeta)*obj.I_0_0*J_omega_0]*zeta ...
                + Ia_0_0*[J_a_zeta_0*zeta-af_0-obj.S(J_omega_0*zeta)*(J_v_0*zeta-vf_0); J_alpha_zeta_0*zeta-obj.S(J_omega_0*zeta)*J_omega_0*zeta] ...
                + [obj.S(J_omega_0*zeta) zeros(3,3); obj.S(J_v_0*zeta-vf_0) obj.S(J_omega_0*zeta)]*Ia_0_0*[J_v_0*zeta-vf_0; J_omega_0*zeta]);
            C1 = [J_v_1; J_omega_1]'*([obj.m1*J_a_zeta_1; I_1_0*J_alpha_zeta_1+obj.S(J_omega_1*zeta)*I_1_0*J_omega_1]*zeta ...
                + Ia_1_0*[J_a_zeta_1*zeta-af_0-obj.S(J_omega_0*zeta)*(J_v_1*zeta-vf_0); J_alpha_zeta_1*zeta-obj.S(J_omega_0*zeta)*J_omega_1*zeta] ...
                + [obj.S(J_omega_1*zeta) zeros(3,3); obj.S(J_v_1*zeta-vf_0) obj.S(J_omega_1*zeta)]*Ia_1_0*[J_v_1*zeta-vf_0; J_omega_1*zeta]);
            C2 = [J_v_2; J_omega_2]'*([obj.m2*J_a_zeta_2; I_2_0*J_alpha_zeta_2+obj.S(J_omega_2*zeta)*I_2_0*J_omega_2]*zeta ...
                + Ia_2_0*[J_a_zeta_2*zeta-af_0-obj.S(J_omega_0*zeta)*(J_v_2*zeta-vf_0); J_alpha_zeta_2*zeta-obj.S(J_omega_0*zeta)*J_omega_2*zeta] ...
                + [obj.S(J_omega_2*zeta) zeros(3,3); obj.S(J_v_2*zeta-vf_0) obj.S(J_omega_2*zeta)]*Ia_2_0*[J_v_2*zeta-vf_0; J_omega_2*zeta]);
            C3 = [J_v_3; J_omega_3]'*([obj.m3*J_a_zeta_3; I_3_0*J_alpha_zeta_3+obj.S(J_omega_3*zeta)*I_3_0*J_omega_3]*zeta ...
                + Ia_3_0*[J_a_zeta_3*zeta-af_0-obj.S(J_omega_0*zeta)*(J_v_3*zeta-vf_0); J_alpha_zeta_3*zeta-obj.S(J_omega_0*zeta)*J_omega_3*zeta] ...
                + [obj.S(J_omega_3*zeta) zeros(3,3); obj.S(J_v_3*zeta-vf_0) obj.S(J_omega_3*zeta)]*Ia_3_0*[J_v_3*zeta-vf_0; J_omega_3*zeta]);

            C = C0 + C1 + C2 + C3;

            % vector of external forces (- Fg - Fb - Ff - Fd)
            vr_1_0_normal = J_v_1*zeta - vf_0 - dot(J_v_1*zeta-vf_0,R_1_0*obj.uvms_kinematics.z_1_1)*R_1_0*obj.uvms_kinematics.z_1_1;
            vr_2_0_normal = J_v_2*zeta - vf_0 - dot(J_v_2*zeta-vf_0,R_2_0*obj.uvms_kinematics.z_2_2)*R_2_0*obj.uvms_kinematics.z_2_2;
            vr_3_0_normal = J_v_3*zeta - vf_0 - dot(J_v_3*zeta-vf_0,R_3_0*obj.uvms_kinematics.z_3_3)*R_3_0*obj.uvms_kinematics.z_3_3;

            Fe0 = [J_v_0; J_omega_0]'*([-obj.m0*g_0+obj.rho*V0*(g_0-af_0); zeros(3,1)] ...
                + [pi/2*obj.rho*obj.Cd0*r0^2*norm(J_v_0*zeta-vf_0)*(J_v_0*zeta-vf_0); zeros(3,1)]);
            Fe1 = [J_v_1; J_omega_1]'*([-obj.m1*g_0+obj.rho*V1*(g_0-af_0); zeros(3,1)] ...
                + obj.rho*obj.Cd1*obj.uvms_kinematics.r1*norm(vr_1_0_normal)*...
                [obj.uvms_kinematics.L1*vr_1_0_normal; 1/2*obj.uvms_kinematics.L1^2*obj.S(R_1_0*obj.uvms_kinematics.z_1_1)*vr_1_0_normal]);
            Fe2 = [J_v_2; J_omega_2]'*([-obj.m2*g_0+obj.rho*V2*(g_0-af_0); zeros(3,1)] ...
                + obj.rho*obj.Cd2*obj.uvms_kinematics.r2*norm(vr_2_0_normal)*...
                [obj.uvms_kinematics.L2*vr_2_0_normal; 1/2*obj.uvms_kinematics.L2^2*obj.S(R_2_0*obj.uvms_kinematics.z_2_2)*vr_2_0_normal]);
            Fe3 = [J_v_3; J_omega_3]'*([-obj.m3*g_0+obj.rho*V3*(g_0-af_0); zeros(3,1)] ...
                + obj.rho*obj.Cd3*obj.uvms_kinematics.r3*norm(vr_3_0_normal)*...
                [obj.uvms_kinematics.L3*vr_3_0_normal; 1/2*obj.uvms_kinematics.L3^2*obj.S(R_3_0*obj.uvms_kinematics.z_3_3)*vr_3_0_normal]);
            
            Fe = Fe0 + Fe1 + Fe2 + Fe3;
        end
        
        %% Volume of Link (assume neutrally buoyant m = rho*V)
        function [V0, V1, V2, V3] = GetVolume(obj)
%             V0 = 4/3*pi*(obj.uvms_kinematics.Lx/2)*(obj.uvms_kinematics.Ly/2)*(obj.uvms_kinematics.Lz/2);
%             V1 = pi*obj.uvms_kinematics.r1^2*obj.uvms_kinematics.L1;
%             V2 = pi*obj.uvms_kinematics.r2^2*obj.uvms_kinematics.L2;
%             V3 = pi*obj.uvms_kinematics.r3^2*obj.uvms_kinematics.L3;
            
            V0 = obj.m0/obj.rho;
            V1 = obj.m1/obj.rho;
            V2 = obj.m2/obj.rho;
            V3 = obj.m3/obj.rho;
        end
        
        %% Inertia Matrix of Link        
        % inertia matrix of link 1 2 3 expressed in frame 0
        function [I_1_0, I_2_0, I_3_0] = GetInertiaMatrix(obj)
            R_1_0 = obj.uvms_kinematics.get_R_1_0();
            R_2_0 = obj.uvms_kinematics.get_R_2_0();
            R_3_0 = obj.uvms_kinematics.get_R_3_0();
            
            I_1_0 = R_1_0*obj.I_1_1*R_1_0';
            I_2_0 = R_2_0*obj.I_2_2*R_2_0';
            I_3_0 = R_3_0*obj.I_3_3*R_3_0';
        end
        
        %% Added Mass Matrix
        % added mass matrix of link 0 expressed in frame 0
        function Ia_0_0 = get_Ia_0_0(obj)          
            % semi-axis of vehicle (a > b)
            a = obj.uvms_kinematics.Lx/2;
            b = (obj.uvms_kinematics.Ly*obj.uvms_kinematics.Lz/4)^(1/2);
            
            e = 1 - (b/a)^2;
            m = 4/3*pi*obj.rho*a*b^2;
            alpha0 = 2*(1-e^2)/e^3 * (1/2*log((1+e)/(1-e)) - e);
            beta0  = 1/e^2 - (1-e^2)/(2*e^3) * log((1+e)/(1-e));
            
            Ia_0_0 = zeros(6,6);
            Ia_0_0(1,1) = -alpha0/(2-alpha0)*m;
            Ia_0_0(2,2) = -beta0/(2-beta0)*m;
            Ia_0_0(3,3) = -beta0/(2-beta0)*m;
            Ia_0_0(4,4) = 0;
            Ia_0_0(5,5) = -1/5 * (b^2-a^2)^2*(alpha0-beta0) / (2*(b^2-a^2)+(b^2+a^2)*(beta0-alpha0)) * m;
            Ia_0_0(6,6) = -1/5 * (b^2-a^2)^2*(alpha0-beta0) / (2*(b^2-a^2)+(b^2+a^2)*(beta0-alpha0)) * m;
        end
        
        % added mass matrix of link 1 expressed in frame 1
        function Ia_1_1 = get_Ia_1_1(obj)
            k1 = obj.rho*pi*obj.uvms_kinematics.r1^2*obj.uvms_kinematics.L1/4;           
            Ia_1_1 = diag([k1 k1 0 k1*obj.uvms_kinematics.L1^2/3 k1*obj.uvms_kinematics.L1^2/3 0]);
        end
        
        % added mass matrix of link 2 expressed in frame 2
        function Ia_2_2 = get_Ia_2_2(obj)
            k2 = obj.rho*pi*obj.uvms_kinematics.r2^2*obj.uvms_kinematics.L2/4;           
            Ia_2_2 = diag([k2 k2 0 k2*obj.uvms_kinematics.L2^2/3 k2*obj.uvms_kinematics.L2^2/3 0]);
        end
        
        % added mass matrix of link 3 expressed in frame 3
        function Ia_3_3 = get_Ia_3_3(obj)
            k3 = obj.rho*pi*obj.uvms_kinematics.r3^2*obj.uvms_kinematics.L3/4;            
            Ia_3_3 = diag([k3 k3 0 k3*obj.uvms_kinematics.L3^2/3 k3*obj.uvms_kinematics.L3^2/3 0]);
        end
        
        % added mass matrix of link 0 1 2 3 expressed in frame 0
        function [Ia_0_0, Ia_1_0, Ia_2_0, Ia_3_0] = GetAddedMassMatrix(obj)
            R_1_0 = obj.uvms_kinematics.get_R_1_0();
            R_2_0 = obj.uvms_kinematics.get_R_2_0();
            R_3_0 = obj.uvms_kinematics.get_R_3_0();
            
            Ia_0_0 = obj.get_Ia_0_0();
            Ia_1_1 = obj.get_Ia_1_1();
            Ia_2_2 = obj.get_Ia_2_2();
            Ia_3_3 = obj.get_Ia_3_3();
            
            Ia_1_0 = zeros(6,6);
            Ia_1_0(1:3,1:3) = R_1_0*Ia_1_1(1:3,1:3)*R_1_0';
            Ia_1_0(4:6,4:6) = R_1_0*Ia_1_1(4:6,4:6)*R_1_0';
            
            Ia_2_0 = zeros(6,6);
            Ia_2_0(1:3,1:3) = R_2_0*Ia_2_2(1:3,1:3)*R_2_0';
            Ia_2_0(4:6,4:6) = R_2_0*Ia_2_2(4:6,4:6)*R_2_0';
            
            Ia_3_0 = zeros(6,6);
            Ia_3_0(1:3,1:3) = R_3_0*Ia_3_3(1:3,1:3)*R_3_0';
            Ia_3_0(4:6,4:6) = R_3_0*Ia_3_3(4:6,4:6)*R_3_0';
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