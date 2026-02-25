%% ================================================
% RK4 Orbit + 3-Axis Reaction Wheel Control
% ================================================

clear all; clc; close all;

%% 1) CONSTANTS
mu = 398600e3;           
Re = 6378e3;             
dt = 1;              
tf = 86400;           
N  = floor(tf/dt);

%% 2) INITIAL ORBIT SETUP
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

%% 3) ATTITUDE SETUP
q = zeros(4,N+1);
q(:,1) = [0;0;0;1];      % initial quaternion

theta = zeros(3,N+1);        % body rates
theta(:,1) = [0;0;0];


w = zeros(3,N+1);        % body rates
w(:,1) = [1e-12;1e-12;1e-12];

wRW = zeros(3,N+1);        % body rates
wRW(:,1) = [1e-12;1e-12;1e-12];

wdot = zeros(3,N+1);        % body rates
wdot(:,1) = [0;0;0];

J_sat = [0.00225049 0.00830119 0.00830119];
J_rw  = [0.00225049e-3 0.00830119e-3 0.00830119e-3];
referenceDirections = [1 0 0; 0 1 0; 0 0 1];

%% 4) CONTROL SETUP
theta_des = [pi/2; pi/4; pi/4];   % desired Euler angles
alpha_rw_max = 1e-4;                % wheel accel limit

%% 5) CUBE GEOMETRY
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

%% 6) FIGURE 1 – INERTIAL FRAME
figure(1); clf; hold on; axis equal; grid on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Inertial Frame')

[Xe,Ye,Ze] = sphere(60);
surf(Re*Xe/1000,Re*Ye/1000,Re*Ze/1000,'FaceColor',[0.3 0.6 1],'EdgeColor','none');
lighting gouraud; camlight headlight

hCube1 = patch('Vertices',V0+r(:,1).','Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);
xlim([-1e5 1e5]); ylim([-1e5 1e5]); zlim([-1e5 1e5]);

%% 7) FIGURE 2 – BODY FRAME
figure(2); clf; hold on; axis equal; grid on;
xlabel('X_b'); ylabel('Y_b'); zlabel('Z_b');
title('Body Frame')

lim = max([Lx Ly Lz]);
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

hCube2 = patch('Vertices',V0,'Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

view(45,30)

%% 8) STORAGE
timestamps = zeros(1,N);

%% DYNAMICS FUNCTIONS
f_orbit = @(r_current,v_current)[ ...
    v_current;
    -mu*r_current/norm(r_current)^3 ];
f_angacc  = @(w_current,J) (-inv(J))*cross(w_current,J*w_current);

%% 9) MAIN LOOP
for k = 2:N

    %% ---------- ORBIT RK4 ----------
    k1 = f_orbit(r(:,k),v(:,k));
    k2 = f_orbit(r(:,k)+0.5*dt*k1(1:3), v(:,k)+0.5*dt*k1(4:6));
    k3 = f_orbit(r(:,k)+0.5*dt*k2(1:3), v(:,k)+0.5*dt*k2(4:6));
    k4 = f_orbit(r(:,k)+dt*k3(1:3), v(:,k)+dt*k3(4:6));

    r(:,k+1) = r(:,k) + (dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v(:,k+1) = v(:,k) + (dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));

    %% ---------- ATTITUDE DYNAMICS ----------
    
    controlTorque_RW = axisWheelTorqueOutput(theta,theta_des,J_sat,dt,1,k);
    L_ext = [0 0 0];
    u = (-controlTorque_RW ./ J_rw);
    k1 = satAttitudeControl( w(:,k), wRW(:,k), u, diag(J_sat), diag(J_rw), referenceDirections, controlTorque_RW, L_ext);
    k2 = satAttitudeControl( w(:,k)+0.5*dt*k1, wRW(:,k)+0.5*dt*u', u, diag(J_sat), diag(J_rw), referenceDirections, controlTorque_RW, L_ext);
    k3 = satAttitudeControl( w(:,k)+0.5*dt*k2, wRW(:,k)+0.5*dt*u', u, diag(J_sat), diag(J_rw), referenceDirections, controlTorque_RW, L_ext);
    k4 = satAttitudeControl( w(:,k)+dt*k3,     wRW(:,k)+dt*u',     u, diag(J_sat), diag(J_rw), referenceDirections, controlTorque_RW, L_ext);
    
    w(:,k+1) = w(:,k) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    wRW(:,k+1) = wRW(:,k) + dt*u';
    
    % Integrate angle 
    theta(:,k+1) = theta(:,k) + w(:,k+1)*dt;
    
    
    %% ---------- CONVERT TO EULER ----------
    angles = theta(:,k+1); 
    
    tol = deg2rad(0.1);
    tol_omega = deg2rad(0.1);

    %% ---------- Initialize new Angle Variables ----------

    angles_new_x = angles(1);
    angles_new_y = angles(2);
    angles_new_z = angles(3);

    w_new_x = w(1,k+1);
    w_new_y = w(2,k+1);
    w_new_z = w(3,k+1);

    
    %% ---------- UPDATE STATE ----------
    w(:,k+1) = [w_new_x w_new_y w_new_z];

    % Rebuild quaternion from updated Euler angles
    Rb = eul2rotm([angles_new_x; angles_new_y; angles_new_z]',"ZYX");
    
    q(:,k+1) = rotm2quat(Rb).';
    
    timestamps(k) = k*dt;

    % Inertial display
    set(hCube1,'Vertices',(Rb*V0.').'*5 + (r(:,k+1)/5).');

    % Body display
    set(hCube2,'Vertices',(Rb*V0.').');

    drawnow limitrate


end
