%% ================================================
% RK4 Orbit + 3-Axis Reaction Wheel Control
% Dual Display (Inertial + Body)
% ================================================

clear all; clc; close all;

%% =================================================
%% 1) CONSTANTS
%% =================================================
mu = 398600e3;           
Re = 6378e3;             
dt = 1;              
tf = 86400;           
N  = floor(tf/dt);

%% =================================================
%% 2) INITIAL ORBIT SETUP
%% =================================================
a  = 400e3;
e  = 0.00007;
inc  = deg2rad(53);
RAAN = deg2rad(219.77);
argp = deg2rad(332.17);
E0   = deg2rad(44.66);

n = sqrt(mu/a^3);

R3 = @(th)[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
R1 = @(th)[1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
A = R3(RAAN)*R1(inc)*R3(argp);

r_pf = a*[cos(E0)-e; sqrt(1-e^2)*sin(E0); 0];
v_pf = (n*a^2/norm(r_pf))*[-sin(E0); sqrt(1-e^2)*cos(E0); 0];

r = zeros(3,N+1);
v = zeros(3,N+1);

r(:,1) = A.'*r_pf;
v(:,1) = A.'*v_pf;

%% =================================================
%% 3) ATTITUDE SETUP
%% =================================================
q = zeros(4,N+1);
q(:,1) = [0;0;0;1];      % initial quaternion

w = zeros(3,N+1);        % body rates
w(:,1) = [1e-12;1e-12;1e-12];

J_sat = diag([0.00225049 0.00830119 0.00830119]);
J_rw  = diag([0.00225049e-3 0.00830119e-3 0.00830119e-3]);

%% =================================================
%% 4) CONTROL SETUP
%% =================================================
theta_des = [pi/3; pi/2; pi/4];   % desired Euler angles
alpha_rw_max = 1e-4;                % wheel accel limit

%% =====================================================
%% 5) CUBE GEOMETRY
%% =====================================================
Lx=400; Ly=400; Lz=800;

V0 = [
   -Lx/2 -Ly/2 -Lz/2;
    Lx/2 -Ly/2 -Lz/2;
    Lx/2  Ly/2 -Lz/2;
   -Lx/2  Ly/2 -Lz/2;
   -Lx/2 -Ly/2  Lz/2;
    Lx/2 -Ly/2  Lz/2;
    Lx/2  Ly/2  Lz/2;
   -Lx/2  Ly/2  Lz/2 ];

F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

%% =====================================================
%% 6) FIGURE 1 – INERTIAL FRAME
%% =====================================================
figure(1); clf; hold on; axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Inertial Frame')

[Xe,Ye,Ze] = sphere(60);
surf(Re*Xe/1000,Re*Ye/1000,Re*Ze/1000,'FaceColor',[0.3 0.6 1],'EdgeColor','none');
lighting gouraud; camlight headlight

hCube1 = patch('Vertices',V0+r(:,1).','Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);
xlim([-1e5 1e5]); ylim([-1e5 1e5]); zlim([-1e5 1e5]);

%% =====================================================
%% 7) FIGURE 2 – BODY FRAME
%% =====================================================
figure(2); clf; hold on; axis equal; grid on;
xlabel('X_b'); ylabel('Y_b'); zlabel('Z_b');
title('Body Frame')

lim = max([Lx Ly Lz]);
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

hCube2 = patch('Vertices',V0,'Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

view(45,30)


%% =================================================
%% 8) STORAGE
%% =================================================
timestamps = zeros(1,N);

%% =================================================
%% DYNAMICS FUNCTIONS
%% =================================================
f_orbit = @(r_current,v_current)[ ...
    v_current;
    -mu*r_current/norm(r_current)^3 ];
f_angacc  = @(w_current,J) (-inv(J))*cross(w_current,J*w_current);

%% =================================================
%% 9) MAIN LOOP
%% =================================================
for k = 1:N

    %% ---------- ORBIT RK4 ----------
    k1 = f_orbit(r(:,k),v(:,k));
    k2 = f_orbit(r(:,k)+0.5*dt*k1(1:3), v(:,k)+0.5*dt*k1(4:6));
    k3 = f_orbit(r(:,k)+0.5*dt*k2(1:3), v(:,k)+0.5*dt*k2(4:6));
    k4 = f_orbit(r(:,k)+dt*k3(1:3), v(:,k)+dt*k3(4:6));

    r(:,k+1) = r(:,k) + (dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v(:,k+1) = v(:,k) + (dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));


    %% ---------- ATTITUDE DYNAMICS ----------
    k1 = f_angacc(w(:,k),J_sat);
    k2 = f_angacc(w(:,k)+0.5*dt*k1,J_sat);
    k3 = f_angacc(w(:,k)+0.5*dt*k2,J_sat);
    k4 = f_angacc(w(:,k)+dt*k3,J_sat);

    w(:,k+1) = w(:,k) + (dt/6)*(k1+2*k2+2*k3+k4);

    %% ---------- CONVERT TO EULER ----------
    angles = quat2eul(q(:,k).');  
    angles = angles(:);
    tol = deg2rad(1);
    tol_omega = deg2rad(0.1);

    %% ---------- Initialize new Angle Variables ----------

    angles_new_x = mod(angles(1),pi);
    angles_new_y = mod(angles(2),pi);
    angles_new_z = mod(angles(3),pi);

    w_new_x = w(1,k+1);
    w_new_y = w(2,k+1);
    w_new_z = w(3,k+1);

    %% ---------- 1-AXIS REACTION WHEEL CONTROLS, Deine new Angles ----------
    if (norm(angles(3)-theta_des(3),2) > tol) % Z axis control
        [angles_new_z, w_new_z] = axisWheelSim( ...
            theta_des(3), ...
            angles(3), ...
            w(3,k+1), ...
            alpha_rw_max, ...
            dt, ...
            J_sat(3,3), ...
            J_rw(3,3),tol_omega);
    elseif (norm(abs(angles(1:2)-theta_des(1:2)),2) > tol) % XY axis control
        [angles_new_x, w_new_x] = axisWheelSim( ...
            theta_des(1), ...
            angles(1), ...
            w(1,k+1), ...
            alpha_rw_max, ...
            dt, ...
            J_sat(1,1), ...
            J_rw(1,1),tol_omega);
    
        [angles_new_y, w_new_y] = axisWheelSim( ...
            theta_des(2), ...
            angles(2), ...
            w(2,k+1), ...
            alpha_rw_max, ...
            dt, ...
            J_sat(2,2), ...
            J_rw(2,2),tol_omega);
    end

    %% ---------- UPDATE STATE ----------
    w(:,k+1) = [w_new_x w_new_y w_new_z];
    disp([abs(([angles_new_x,angles_new_y,angles_new_z]'-theta_des)), [angles_new_x,angles_new_y,angles_new_z]', theta_des]);

    % Rebuild quaternion from updated Euler angles
    Rb = eul2rotm([angles_new_x; angles_new_y; angles_new_z].',"ZYX");
    q(:,k+1) = rotm2quat(Rb).';
    
    timestamps(k) = k*dt;

    % Inertial display
    set(hCube1,'Vertices',(Rb*V0.').'*5 + (r(:,k+1)/5).');

    % Body display
    set(hCube2,'Vertices',(Rb*V0.').');

    drawnow limitrate


end
