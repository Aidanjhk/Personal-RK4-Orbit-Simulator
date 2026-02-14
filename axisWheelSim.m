function [theta, omega, thetaOld] = axisWheelSim( ...
    theta, omega, thetaOld, ...
    theta_des, dt, wn, zeta, I, Td)
    % INPUT:
    % theta --> current rotation vector
    % omega --> current angular velocity vector
    % thetaOld --> previous rotation vector
    % theta_des --> desired rotation vector
    % dt --> time step
    % wn --> natural undamped frequency
    % zeta --> damping ratio
    % I --> inertia diagonal  
    % Td --> external & disturbance torque
    % 
    % OUTPUT:
    % theta --> new rotation vector
    % omega --> new angular velocity vector
    % thetaOld --> current rotation vector

    Ivec = I(:);
    
    Kf = (wn^2).*Ivec;
    tau_r = (2*zeta*wn.*Ivec)./Kf;
    
    theta_e = theta - theta_des;
    
    Tc = -Kf .* ( theta_e + tau_r .* ((theta_e - thetaOld)/dt) );
    
    omegaDot = (Tc + Td)./Ivec;
    omega = omega + dt.*omegaDot;
    theta = theta + dt.*omega;
    
    thetaOld = theta_e;
end