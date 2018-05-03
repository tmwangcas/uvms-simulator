classdef UvmsGraphics
    
    properties
        uvms_dynamics;
    end
    
    methods
        function obj = UvmsGraphics(uvms_dynamics)
            if nargin == 1
                obj.uvms_dynamics = uvms_dynamics;
            else
                error('Not enough parameters');
            end
        end
        
        %% Update UVMS Dynamics
        function obj = UpdateUvmsDynamics(obj, uvms_dynamics)
            obj.uvms_dynamics = uvms_dynamics;
        end
        
        %% Draw Robot
        % draw vehicle
        function DrawVehicle(obj)
            vehicle_vertex = obj.uvms_dynamics.uvms_kinematics.GetVehicleGeometry();
            
            % color definition
            orange = [1    0.5  0];
            gray   = [0.6  0.6  0.6];
            blue   = [0.25 0.4  0.9];
            
            % draw up & down surfaces
            ff1 = fill3(vehicle_vertex(1,[1 2 10 9 1]),vehicle_vertex(2,[1 2 10 9 1]),vehicle_vertex(3,[1 2 10 9 1]),orange);
            ff2 = fill3(vehicle_vertex(1,[3 4 12 11 3]),vehicle_vertex(2,[3 4 12 11 3]),vehicle_vertex(3,[3 4 12 11 3]),orange);
            ff3 = fill3(vehicle_vertex(1,[13 14 15 16 13]),vehicle_vertex(2,[13 14 15 16 13]),vehicle_vertex(3,[13 14 15 16 13]),orange);
            ff4 = fill3(vehicle_vertex(1,[5 6 7 8 5]),vehicle_vertex(2,[5 6 7 8 5]),vehicle_vertex(3,[5 6 7 8 5]),gray);
            
            % draw front & back surfaces
            ff5 = fill3(vehicle_vertex(1,[1 5 8 4 12 16 13 9 1]),vehicle_vertex(2,[1 5 8 4 12 16 13 9 1]),vehicle_vertex(3,[1 5 8 4 12 16 13 9 1]),blue);
            ff6 = fill3(vehicle_vertex(1,[2 6 7 3 11 15 14 10 2]),vehicle_vertex(2,[2 6 7 3 11 15 14 10 2]),vehicle_vertex(3,[2 6 7 3 11 15 14 10 2]),blue);
            
            % draw left & right surfaces
            ff7 = fill3(vehicle_vertex(1,[1 2 6 5 1]),vehicle_vertex(2,[1 2 6 5 1]),vehicle_vertex(3,[1 2 6 5 1]),gray);        
            ff8 = fill3(vehicle_vertex(1,[3 4 8 7 3]),vehicle_vertex(2,[3 4 8 7 3]),vehicle_vertex(3,[3 4 8 7 3]),gray);        
            ff9 = fill3(vehicle_vertex(1,[9 10 14 13 9]),vehicle_vertex(2,[9 10 14 13 9]),vehicle_vertex(3,[9 10 14 13 9]),gray);        
            ff10 = fill3(vehicle_vertex(1,[12 11 15 16 12]),vehicle_vertex(2,[12 11 15 16 12]),vehicle_vertex(3,[12 11 15 16 12]),gray); 
            
            % set transparency
%             set(ff1,'FaceAlpha',1);
%             set(ff2,'FaceAlpha',1);
%             set(ff3,'FaceAlpha',1);
%             set(ff4,'FaceAlpha',1);
%             set(ff5,'FaceAlpha',1);
%             set(ff6,'FaceAlpha',1);
%             set(ff7,'FaceAlpha',1);
%             set(ff8,'FaceAlpha',1);
%             set(ff9,'FaceAlpha',1);
%             set(ff10,'FaceAlpha',1);
        end
        
        % draw manipulator
        function DrawManipulator(obj)
            [po_1_I, po_2_I, po_3_I, p_L2_left_I, p_L2_right_I, p_L2_start_I, p_L2_end_I] = ...
                obj.uvms_dynamics.uvms_kinematics.GetManipulatorGeometry();
            
            % draw link 1
            obj.DrawLink(po_1_I, po_2_I, obj.uvms_dynamics.uvms_kinematics.r1);
            obj.DrawLink(p_L2_left_I, p_L2_right_I, obj.uvms_dynamics.uvms_kinematics.r1);
            
            % draw link 2
            obj.DrawLink(p_L2_start_I, p_L2_end_I, obj.uvms_dynamics.uvms_kinematics.r2);
            
            % draw link 3
            obj.DrawLink(p_L2_end_I, po_3_I, obj.uvms_dynamics.uvms_kinematics.r3);
        end
        
        % draw desired trajectory of end-effector
        function DrawEndEffectorDesiredTrajectory(obj)            
            cx = obj.uvms_dynamics.uvms_kinematics.p_ee_init(1);
            cy = obj.uvms_dynamics.uvms_kinematics.p_ee_init(2);
            cz = obj.uvms_dynamics.uvms_kinematics.p_ee_init(3) + obj.uvms_dynamics.uvms_kinematics.r_ee_traj;
            
            alpha = 0:pi/1000:2*pi;
            
            x = 0*alpha + cx;
            y = obj.uvms_dynamics.uvms_kinematics.r_ee_traj*cos(alpha) + cy;
            z = obj.uvms_dynamics.uvms_kinematics.r_ee_traj*sin(alpha) + cz;
            
            plot3(x,y,z,'Color',[0 0.8 0],'LineWidth',3);
        end
        
        %% Plot Variables
        % draw end-effector trajectory
        function DrawEndEffectorTracjectory(obj, p_ee, p_ee_d)
            figure
            set(gcf,'outerposition',get(0,'screensize'))

            plot3(p_ee_d(1,:),p_ee_d(2,:),p_ee_d(3,:),'b','LineWidth',2)
            hold on
            plot3(p_ee(1,:),p_ee(2,:),p_ee(3,:),'Color',[0 0.8 0],'LineWidth',4);
            hold on

            plot3(p_ee(1,1),p_ee(2,1),p_ee(3,1),'ok','MarkerFaceColor','y','MarkerSize',10)
            hold on
            plot3(p_ee(1,end),p_ee(2,end),p_ee(3,end),'sk','MarkerFaceColor','r','MarkerSize',10)

            grid on
            title('3D Trajectories')
            legend('desired','actual')
            xlabel('X(m)')
            ylabel('Y(m)')
            zlabel('Z(m)')
            view(90,0)
            axis([0.276-0.3 0.276+0.3 -0.3 0.3 0.37+0.24-0.3 0.37+0.24+0.3])
        end
        
        % plot end-effector trajectory
        function PlotEndEffectorTracjectory(obj, t_vec, p_ee, p_ee_d)
            figure
            set(gcf,'outerposition',get(0,'screensize'))

            % p_ee_d
            subplot(3,2,1)
            plot(t_vec,p_ee_d(1,:))
            grid on
            title('p\_ee\_d')
            xlabel('time/s')
            subplot(3,2,3)
            plot(t_vec,p_ee_d(2,:))
            grid on
            xlabel('time/s')
            subplot(3,2,5)
            plot(t_vec,p_ee_d(3,:))
            grid on
            xlabel('time/s')

            % p_ee
            subplot(3,2,2)
            plot(t_vec,p_ee(1,:))
            grid on
            title('p\_ee')
            xlabel('time/s')
            subplot(3,2,4)
            plot(t_vec,p_ee(2,:))
            grid on
            xlabel('time/s')
            subplot(3,2,6)
            plot(t_vec,p_ee(3,:))
            grid on
            xlabel('time/s')
        end
        
        % plot generalized coordinate
        function PlotGeneralizedCoordinate(obj, t_vec, eta)
            figure
            set(gcf,'outerposition',get(0,'screensize'))

            % p
            subplot(3,3,1)
            plot(t_vec,eta(1,:))
            grid on
            title('p')
            xlabel('time/s')
            subplot(3,3,4)
            plot(t_vec,eta(2,:))
            grid on
            xlabel('time/s')
            subplot(3,3,7)
            plot(t_vec,eta(3,:))
            grid on
            xlabel('time/s')

            % theta
            subplot(3,3,2)
            plot(t_vec,eta(4,:))
            grid on
            title('theta')
            xlabel('time/s')
            subplot(3,3,5)
            plot(t_vec,eta(5,:))
            grid on
            xlabel('time/s')
            subplot(3,3,8)
            plot(t_vec,eta(6,:))
            grid on
            xlabel('time/s')

            % q
            subplot(3,3,3)
            plot(t_vec,eta(7,:))
            grid on
            title('q')
            xlabel('time/s')
            subplot(3,3,6)
            plot(t_vec,eta(8,:))
            grid on
            xlabel('time/s')
            subplot(3,3,9)
            plot(t_vec,eta(9,:))
            grid on
            xlabel('time/s')           
        end
        
        % plot generalized velocity
        function PlotGeneralizedVelocity(obj, t_vec, zeta)
            figure
            set(gcf,'outerposition',get(0,'screensize'))

            % v
            subplot(3,3,1)
            plot(t_vec,zeta(1,:))
            grid on
            title('v')
            xlabel('time/s')
            subplot(3,3,4)
            plot(t_vec,zeta(2,:))
            grid on
            xlabel('time/s')
            subplot(3,3,7)
            plot(t_vec,zeta(3,:))
            grid on
            xlabel('time/s')

            % omega
            subplot(3,3,2)
            plot(t_vec,zeta(4,:))
            grid on
            title('omega')
            xlabel('time/s')
            subplot(3,3,5)
            plot(t_vec,zeta(5,:))
            grid on
            xlabel('time/s')
            subplot(3,3,8)
            plot(t_vec,zeta(6,:))
            grid on
            xlabel('time/s')

            % dq
            subplot(3,3,3)
            plot(t_vec,zeta(7,:))
            grid on
            title('dq')
            xlabel('time/s')
            subplot(3,3,6)
            plot(t_vec,zeta(8,:))
            grid on
            xlabel('time/s')
            subplot(3,3,9)
            plot(t_vec,zeta(9,:))
            grid on
            xlabel('time/s')           
        end
        
        % plot control forces
        function PlotControlForces(obj, t_vec, tau_c)
            figure
            set(gcf,'outerposition',get(0,'screensize'))

            % vehicle force
            subplot(3,3,1)
            plot(t_vec,tau_c(1,:))
            grid on
            title('vehicle force')
            xlabel('time/s')
            subplot(3,3,4)
            plot(t_vec,tau_c(2,:))
            grid on
            xlabel('time/s')
            subplot(3,3,7)
            plot(t_vec,tau_c(3,:))
            grid on
            xlabel('time/s')

            % vehicle torque
            subplot(3,3,2)
            plot(t_vec,tau_c(4,:))
            grid on
            title('vehicle torque')
            xlabel('time/s')
            subplot(3,3,5)
            plot(t_vec,tau_c(5,:))
            grid on
            xlabel('time/s')
            subplot(3,3,8)
            plot(t_vec,tau_c(6,:))
            grid on
            xlabel('time/s')

            % joint torque
            subplot(3,3,3)
            plot(t_vec,tau_c(7,:))
            grid on
            title('joint torque')
            xlabel('time/s')
            subplot(3,3,6)
            plot(t_vec,tau_c(8,:))
            grid on
            xlabel('time/s')
            subplot(3,3,9)
            plot(t_vec,tau_c(9,:))
            grid on
            xlabel('time/s')           
        end
    end
    
    methods(Static)
        % draw cylinder link
        function DrawLink(pa, pb, r)
            % input:
            %   pa   dim 3x1    "starting" point
            %   pb   dim 3x1    "ending" point
            %   r    dim 1x1    link radius

            % generate points on the cylinder aligned with x
            [z, y, x] = cylinder(r*ones(41,1), 40);
            x = norm(pb - pa)*x;
            x2 = (pb - pa)/norm(pb - pa);
            if ((x2(1) == -1) || (x2(1) == 1))
                z2 = [0; 0; 1];
            else
                z2 = cross([1; 0; 0], x2); 
                z2 = z2/norm(z2);
            end
            y2 = cross(x2, z2);

            % rotate and translate points
            R = [x2 y2 z2];
            for i = 1:length(x)
                for j = 1:length(x)
                    % rotation
                    rr = R*[x(i,j); y(i,j); z(i,j)];
                    x(i,j) = rr(1);
                    y(i,j) = rr(2); 
                    z(i,j) = rr(3);
                    % translation
                    x(i,j) = x(i,j) + pa(1);
                    y(i,j) = y(i,j) + pa(2); 
                    z(i,j) = z(i,j) + pa(3);
                end
            end
            
            h = surfl(x, y, z, [0 0 -5]);
            set(h, 'facealpha', 0.5)
            set(h, 'facecolor', 'interp');
            set(h, 'edgecolor', 'none');
            colormap(bone)
        end
    end
end