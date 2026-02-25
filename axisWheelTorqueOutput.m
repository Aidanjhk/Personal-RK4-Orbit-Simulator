function [Torque] = axisWheelTorqueOutput(anglesCurrent,anglesDesired,SatInertia,dt,zeta,k)
    % INPUT:
    % 
    % 
    % OUTPUT:
    % 
    % 

    Ivec = [SatInertia(1),SatInertia(2),SatInertia(3)];
    wn = 0.01;
    
    Kf = (wn^2)*Ivec;
    tau_r = (2*zeta*wn.*Ivec)./Kf;
    anglesCurrent(:,k-1)
    theta_e_prev = anglesCurrent(:,k-1) - anglesDesired;
    theta_e = anglesCurrent(:,k) - anglesDesired;
    
    Tc = -Kf * ( theta_e + (tau_r * (theta_e - theta_e_prev)/dt) );
    Torque = Tc;
end