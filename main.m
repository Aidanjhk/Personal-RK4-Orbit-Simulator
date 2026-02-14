%% ================================
% RK4 Orbit + Attitude with
% Dual Display (Inertial + Target-Centered)
% ================================

clear; clc; close all;

%% ---------- Constants ----------
mu = 398600;           % km^3/s^2
Re = 6378;             % km
dt = 0.1;                % s
tf = 34*3600;           % 3 hours
N  = floor(tf/dt);

%% ---------- Initial Orbit ----------
a  = 35790;            % km
e  = 0.00007;
inc  = deg2rad(0.04);
RAAN = deg2rad(219.77);
argp = deg2rad(332.17);
E0   = deg2rad(44.66);

n = sqrt(mu/a^3);

% Rotation matrices
R3 = @(th)[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
R1 = @(th)[1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
A = R3(RAAN)*R1(inc)*R3(argp);

% Initial state
r_pf = a*[cos(E0)-e; sqrt(1-e^2)*sin(E0); 0];
v_pf = (n*a^2/norm(r_pf))*[-sin(E0); sqrt(1-e^2)*cos(E0); 0];

r = A.'*r_pf;
v = A.'*v_pf;

%% ---------- Attitude ----------
q = [0;0;0;1];               % quaternion
w = [0.005;0.004;0.004];     % rad/s
J = diag([0.25 0.3 0.2]);

%% ---------- Storage ----------
positions = zeros(3,N);
quats = zeros(4,N);

%% ---------- Dynamics ----------
f_orbit = @(r,v)[v; -mu*r/norm(r)^3];
f_w = @(w) -J\(cross(w,J*w));
f_q = @(q,w) 0.5 * [ ...
     0   -w.'; 
     w   -skew(w)] * q;

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

%% ---------- FIGURE 1: Inertial ----------
figure(1); clf; hold on; axis equal; grid on;
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title('Inertial Frame (Orbit)')

[Xe,Ye,Ze] = sphere(60);
surf(Re*Xe,Re*Ye,Re*Ze,'FaceColor',[0.3 0.6 1],'EdgeColor','none');
lighting gouraud; camlight headlight

hCube1 = patch('Vertices',V0+r.','Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

xlim([-5e4 5e4]); ylim([-5e4 5e4]); zlim([-5e4 5e4]);

%% ---------- FIGURE 2: Target-Centered ----------
figure(2); clf; hold on; axis equal; grid on;
xlabel('X_b'); ylabel('Y_b'); zlabel('Z_b');
title('Target-Centered (Body View)')
lim = max([Lx Ly Lz]);
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);

hCube2 = patch('Vertices',V0,'Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);

view(45,30)

%% ---------- RK4 + ANIMATION ----------
for k = 1:N

    % Store
    positions(:,k) = r;
    quats(:,k) = q;

    % --- RK4 Orbit ---
    k1 = f_orbit(r,v);
    k2 = f_orbit(r+0.5*dt*k1(1:3), v+0.5*dt*k1(4:6));
    k3 = f_orbit(r+0.5*dt*k2(1:3), v+0.5*dt*k2(4:6));
    k4 = f_orbit(r+dt*k3(1:3), v+dt*k3(4:6));

    r = r + dt/6*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v = v + dt/6*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));

    % --- Euler Attitude ---
    w = w + dt*f_w(w);
    q = q + dt*f_q(q,w);
    q = q/norm(q);

    % --- Rotation matrix ---
    Rb = quat2dcm(q.');

    % --- Update displays ---
    set(hCube1,'Vertices',(Rb*V0.').' + r.');
    set(hCube2,'Vertices',(Rb*V0.').');


    drawnow limitrate
end

%% ---------- Utility ----------
function S = skew(w)
S = [  0   -w(3)  w(2);
      w(3)   0   -w(1);
     -w(2)  w(1)   0  ];
end
