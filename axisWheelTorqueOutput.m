function [Torque] = axisWheelTorqueOutput(anglesCurrent,anglesDesired,SatInertia,dt,zeta,k)
    % INPUT:
    % 
    % 
    % OUTPUT:
    % 
    % 

    Ivec = [SatInertia(1),SatInertia(2),SatInertia(3)];
    wn = 0.01;
    
    Kf = (wn^2).*Ivec;
    tau_r = (2*zeta*wn*Ivec)./Kf;
    anglesCurrent(:,k-1)
    theta_e_prev = anglesCurrent(:,k-1) - anglesDesired;
    theta_e = anglesCurrent(:,k) - anglesDesired;
    
    Tc(1) = -Kf(1) * ( theta_e(1) + (tau_r(1) .* (theta_e(1) - theta_e_prev(1))/dt) );
    Tc(2) = -Kf(2) * ( theta_e(2) + (tau_r(2) .* (theta_e(2) - theta_e_prev(2))/dt) );
    Tc(3) = -Kf(3) * ( theta_e(3) + (tau_r(3) .* (theta_e(3) - theta_e_prev(3))/dt) );
    Torque = Tc;
end