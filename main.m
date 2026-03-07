%% ================================================
% RK4 Orbit + 3-Axis Reaction Wheel Control
% ================================================

clear all; clc; close all;

%% 1) CONSTANTS
mu = 3.98600e14; 
Acs = 0.03;
angle = 0;
Re = 6378e3;  
m = 5;
dt = 60;              
tf = 3600*24*365;           
N  = floor(tf/dt);

%% 2) INITIAL ORBIT SETUP
a  = 450e3 + Re;
e  = 0.0005;
inc  = deg2rad(90);
RAAN = deg2rad(90);
argp = deg2rad(0);
E0   = deg2rad(0);

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
theta_des = [pi/2; pi/3; pi/4];   % desired Euler angles
alpha_rw_max = 1e-4;                % wheel accel limit

%% 5) CUBE GEOMETRY
Lx=400; Ly=400; Lz=1200;

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
surf(Re*Xe/10,Re*Ye/10,Re*Ze/10,'FaceColor',[0.3 0.6 1],'EdgeColor','none');
lighting gouraud; camlight headlight

hCube1 = patch('Vertices',V0+r(:,1).','Faces',F,...
    'FaceColor','r','EdgeColor','k','FaceAlpha',0.9);
xlim([-1e7 1e7]); ylim([-1e7 1e7]); zlim([-1e7 1e7]);

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
f_orbit = @(r_current,v_current,adrag)[ ...
    v_current;
    (-mu*r_current/norm(r_current,2)^3) + (adrag)];
f_angacc  = @(w_current,J) (-inv(J))*cross(w_current,J*w_current);

%% 9) MAIN LOOP
for k = 1:N
    [Tatm,~,~,rho] = atmosisa(norm(r(:,k),2));
    norm(r(:,k),2)
    sigma_t = @(theta) 0.93 - (1.48e-3)*theta - (7e-5)*(theta^2);
    sigma_n = @(theta) 0.63*(1 - exp(-3.38e-2 * theta));

    
    %% ---------- ORBIT RK4 ----------
    Cd_k1 = dragCoefficientCalculator(8314.5,Tatm, 26.98, norm(v(:,k),2), sigma_t, sigma_n, deg2rad(angle),r(:,k));
    Fd_k1 = dragCalcRectangularPrism(Acs, deg2rad(angle), rho, Cd_k1, v(:,k));
    a_drag_k1 = Fd_k1/m;
    k1 = f_orbit(r(:,k),v(:,k),a_drag_k1);
    
    Cd_k2 = dragCoefficientCalculator(8314.5,Tatm, 26.98, norm(v(:,k),2), sigma_t, sigma_n, deg2rad(angle),r(:,k)+0.5*dt*k1(1:3));
    Fd_k2 = dragCalcRectangularPrism(Acs, deg2rad(angle), rho, Cd_k2, v(:,k)+0.5*dt*k1(4:6));
    a_drag_k2 = Fd_k2/m;
    k2 = f_orbit(r(:,k)+0.5*dt*k1(1:3), v(:,k)+0.5*dt*k1(4:6),a_drag_k2);

    Cd_k3 = dragCoefficientCalculator(8314.5,Tatm, 26.98, norm(v(:,k),2), sigma_t, sigma_n, deg2rad(angle),r(:,k)+0.5*dt*k2(1:3));
    Fd_k3 = dragCalcRectangularPrism(Acs, deg2rad(angle), rho, Cd_k3, v(:,k)+0.5*dt*k2(4:6));
    a_drag_k3 = Fd_k3/m;
    k3 = f_orbit(r(:,k)+0.5*dt*k2(1:3), v(:,k)+0.5*dt*k2(4:6),a_drag_k3);

    Cd_k4 = dragCoefficientCalculator(8314.5,Tatm, 26.98, norm(v(:,k),2), sigma_t, sigma_n, deg2rad(angle),r(:,k)+dt*k3(1:3));
    Fd_k4 = dragCalcRectangularPrism(Acs, deg2rad(angle), rho, Cd_k4, v(:,k)+dt*k3(4:6));
    a_drag_k4 = Fd_k4/m;
    k4 = f_orbit(r(:,k)+dt*k3(1:3), v(:,k)+dt*k3(4:6),a_drag_k4);

    r(:,k+1) = r(:,k) + (dt/6)*(k1(1:3)+2*k2(1:3)+2*k3(1:3)+k4(1:3));
    v(:,k+1) = v(:,k) + (dt/6)*(k1(4:6)+2*k2(4:6)+2*k3(4:6)+k4(4:6));

    %% ---------- ATTITUDE DYNAMICS ----------
    if (k > 1)
        controlTorque_RW = axisWheelTorqueOutput(theta,theta_des,J_sat,dt,1/sqrt(2),k);
    else
        controlTorque_RW = 0;
    end
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
    set(hCube1,'Vertices',(Rb*V0.').'*5*100 + (r(:,k+1)/5).');

    % Body display
    set(hCube2,'Vertices',(Rb*V0.').');
    
    drawnow
    


end

figure(3)

plot3(r(1,:),r(2,:),r(3,:))
