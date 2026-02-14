%% ================================
% RK4 Orbit + 3-Axis Wheel Control
% Dual Display (Inertial + Body)
% ================================

clear; clc; close all;

%% ---------- Constants ----------
mu = 398600;           
Re = 6378;             
dt = 0.1;              
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

%% ---------- Orbit Dynamics ----------
f_orbit = @(r,v)[v; -mu*r/norm(r)^3];

%% ---------- Cube Geometry ----------
Lx=800; Ly=400; Lz=400;

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

%% ---------- MAIN LOOP ----------
for k = 1:N

    positions(:,k) = r;
    quats(:,k)     = q;

    % ---- Orbit RK4 ----
    k1 = f_orbit(r,v);
    k2 = f_orbit(r+0.5*dt*k1(1:3), v+0.5*dt*k1(4:6));
    k3 = f_orbit(r+0.5*dt*k2(1:3), v+0.5*dt*k2(4:6));
    k4 = f_orbit(r+dt*k3(1:3), v+dt*k3(4:6));

    r = r + (dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v = v + (dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));

    % ---- Rigid Body Dynamics ----
    w = w + dt * (-J\(cross(w,J*w)));

    % ---- Quaternion Increment ----
    angle = norm(w)*dt;

    axis = w / norm(w);
    dq = axang2quat([axis.' angle]);
   

    q_mat = [q(4) q(1) q(2) q(3)];
    q_mat = quatmultiply(q_mat,dq);
    q_mat = q_mat / norm(q_mat);
    q = [q_mat(2); q_mat(3); q_mat(4); q_mat(1)];

    % ---- Extract Rotation Matrix from Quaternions ----
    Rb = quat2dcm([q(4) q(1) q(2) q(3)]);

    phi = atan2(Rb(3,2),Rb(3,3));
    s   = max(-1,min(1,Rb(3,1)));
    theta = -asin(s);
    psi = atan2(Rb(2,1),Rb(1,1));

    angles = [phi;theta;psi];

    % ---- Desired Angular Position ----
    angles_des = [0; pi/4; pi/4];

    % ---- 3-axis controller ----
    [angles,w,thetaOld] = axisWheelSim( ...
        angles,w,thetaOld, ...
        angles_des,dt, ...
        0.50,1.0,diag(J),[1e-4;1e-4;1e-4]);

    x_err = angles_des(1) - angles(1);
    y_err = angles_des(2) - angles(2);
    z_err = angles_des(2) - angles(2);
    disp([x_err y_err]);


    tol = 1e-4;

    % ---- Build rotation matrix from controlled angles ----
    Rb = Rx(angles(1))*Ry(angles(2));
    

    % ---- Update displays ----
    set(hCube1,'Vertices',(Rb*V0.').' + r.');
    set(hCube2,'Vertices',(Rb*V0.').');

    drawnow limitrate
end



%% Rotation helpers
function R = Rx(a)
R = [1 0 0;0 cos(a) -sin(a);0 sin(a) cos(a)];
end

function R = Ry(a)
R = [cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
end

function R = Rz(a)
R = [cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1];
end
