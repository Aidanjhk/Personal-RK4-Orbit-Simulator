%% ================================
% RK4 Orbit + 3-Axis Wheel Control
% Dual Display (Inertial + Body)
% ================================

clear; clc; close all;


f_quat  = @(q,w) 0.5*w*q;
%% ---------- Constants ----------
mu = 398600;           
Re = 6378;             
dt = 0.01;              
tf = 86400;              % shorten for testing
N  = floor(tf/dt);

%% ---------- Initial Orbit ----------
a  = 35790;
e  = 0.00007;
inc  = deg2rad(0.04);
RAAN = deg2rad(219.77);
argp = deg2rad(332.17);
E0   = deg2rad(44.66);

n = sqrt(mu/a^3);

R3 = @(th)[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
R1 = @(th)[1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
A = R3(RAAN)*R1(inc)*R3(argp);

r_pf = a*[cos(E0)-e; sqrt(1-e^2)*sin(E0); 0];
v_pf = (n*a^2/norm(r_pf))*[-sin(E0); sqrt(1-e^2)*cos(E0); 0];

r = A.'*r_pf;
v = A.'*v_pf;

%% ---------- Attitude ----------
q = [0;0;0;1];       
w = [0.005;0.004;0.004];
J = diag([0.25 0.3 0.2]);
J_RW = diag([1e-4 1e-4 1e-4]);

thetaOld = zeros(3,1);

%% ---------- Storage ----------
positions = zeros(3,N);
quats     = zeros(4,N);
w = zeros(3,N+1);

%% ---------- Orbit Dynamics ----------
f_orbit = @(r,v)[v; -mu*r/norm(r)^3];

%% ---------- Cube Geometry ----------
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

%% ---------- FIGURE 1 ----------
figure(1); clf; hold on; axis equal; grid on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Inertial Frame')

[Xe,Ye,Ze] = sphere(60);
surf(Re*Xe,Re*Ye,Re*Ze,'FaceColor',[0.3 0.6 1],'EdgeColor','none');
lighting gouraud; camlight headlight

hCube1 = patch('Vertices',V0+r.','Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

xlim([-5e4 5e4]); ylim([-5e4 5e4]); zlim([-5e4 5e4]);

%% ---------- FIGURE 2 ----------
figure(2); clf; hold on; axis equal; grid on;
xlabel('X_b'); ylabel('Y_b'); zlabel('Z_b');
title('Body Frame')

lim = max([Lx Ly Lz]);
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

hCube2 = patch('Vertices',V0,'Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

view(45,30)

%% Helper Functions
function A = EulerAngles(e1,e2,e3,theta1,theta2,theta3)
    I = eye(3);
    v = @(e)[0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
    An = @(e,theta) I - sin(theta)*v(e) + (1-cos(theta))*(v(e)^2); 

    A = An(e1,theta1)*An(e2,theta2)*An(e3,theta3);
end
function q = quaternionCalculator(A)
    q4 = sqrt((trace(A) + 1)/4);
    q3 = (A(2,3)-A(3,2))/(4*q4);
    q2 = (A(3,1)-A(1,3))/(4*q4);
    q1 = (A(1,2)-A(2,1))/(4*q4);

    q = [q1;q2;q3;q4];
end

f_angacc  = @(w) (-inv(J))*cross(w,J*w);


%% Initial Conditions
w(1:3,1) = [0.01;0.01;0.01];
q(1:4,1) = [0;0;0;1];

q_prev = [0;0;0;1];
%% ---------- MAIN LOOP ----------
for k = 1:N

    positions(:,k) = r;
    q(:,k)     = q_prev;

    % ---- Orbit RK4 ----
    k1 = f_orbit(r,v);
    k2 = f_orbit(r+0.5*dt*k1(1:3), v+0.5*dt*k1(4:6));
    k3 = f_orbit(r+0.5*dt*k2(1:3), v+0.5*dt*k2(4:6));
    k4 = f_orbit(r+dt*k3(1:3), v+dt*k3(4:6));

    r = r + (dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v = v + (dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));
    
    k1 = f_angacc(w(:,k));
    k2 = f_angacc(w(:,k) + 0.5*dt*k1(1:3));
    k3 = f_angacc(w(:,k) + 0.5*dt*k2(1:3));
    k4 = f_angacc(w(:,k) + dt*k3(1:3));
    w(:,k+1) = w(:,k) + (dt/6)*(k1(1:3) + 2*k2(1:3) + 2*k3(1:3) + k4(1:3));

    w_matrix= [0 w(3,k+1) -w(2,k+1) w(1,k+1); -w(3,k+1) 0 w(1,k+1) w(2,k+1); w(2,k+1) -w(1,k+1) 0 w(3,k+1); -w(1,k+1) -w(2,k+1) -w(3,k+1) 0];
    q_matrix = [q(4,k) q(3,k) -q(2,k); -q(3,k) q(4,k) q(1,k); q(2,k) -q(1,k) q(4,k); -q(1,k) -q(2,k) -q(3,k)];

    k1 = f_quat(q_matrix,w_matrix);
    k1_processed = [k1(2,3); k1(3,1); k1(1,2); k1(1,1)];
    k2 = f_quat(q_matrix + 0.5*dt*k1,w_matrix);
    k2_processed = [k2(2,3); k2(3,1); k2(1,2); k2(1,1)];
    k3 = f_quat(q_matrix + 0.5*dt*k2,w_matrix);
    k3_processed = [k3(2,3); k3(3,1); k3(1,2); k3(1,1)];
    k4 = f_quat(q_matrix + dt*k3,w_matrix);
    k4_processed = [k4(2,3); k4(3,1); k4(1,2); k4(1,1)];

    quats(:,k+1) = q(:,k) + (dt/6)*(k1_processed + 2*k2_processed + 2*k3_processed + k4_processed);
    q = quats(:,k+1);
    
    timestamps(k) = dt*k;


    % ---- Extract euler angles from Quaternions ----
   
    w_euler = [w(1,k+1) w(2,k+1) w(3,k+1)];
    
    % ---- Desired Angular Position ----
    angles_des = [pi/4 0 pi/2];
    angles = quat2eul([q(1) q(2) q(3) q(4)]);
    
    
    % ---- 3-axis controller ----
    [angles,w_euler,thetaOld] = axisWheelSim( ...
        angles,w_euler,thetaOld, ...
        angles_des',dt, ...
        1/dt,0.50,diag(J),[1e-4;1e-4;1e-4]);
    angles = [angles(1,1) angles(2,2) angles(3,3)];
    x_err = angles_des(1) - angles(1);
    y_err = angles_des(2) - angles(2);
    z_err = angles_des(3) - angles(3);
    disp([x_err y_err z_err]);

    % ---- Build rotation matrix from controlled angles ----
    Rb = EulerAngles([1 0 0],[0 1 0],[0 0 1],angles(1),angles(2),angles(3));

    % ---- Update displays ----
    set(hCube1,'Vertices',(Rb*V0.').' + r.');
    set(hCube2,'Vertices',(Rb*V0.').');

    q_prev = [q(1) q(2) q(3) q(4)];
   
    drawnow limitrate
end


