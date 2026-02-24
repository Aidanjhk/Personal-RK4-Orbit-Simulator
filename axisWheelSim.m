function [theta_new, omega_new] = axisWheelSim( ...
    theta_des, ...     % desired angle [rad]
    theta, ...         % current angle [rad]
    omega, ...         % current angular velocity [rad/s]
    alpha_max_rw, ...  % max wheel angular acceleration [rad/s^2]
    dt, ...
    I_sat, ...
    I_rw, ...
    omega_max)

    % --------------------------------------------
    % Reaction wheel control with braking logic
    % --------------------------------------------
    
    % Angle error
    err = theta_des - theta;

    % Convert wheel acceleration limit to satellite limit
    alpha_max_sat = abs((I_rw ./ I_sat) * alpha_max_rw);
    omega_max_sat = alpha_max_sat*dt;
    
    % Direction to target
    dir = sign(err);
    if dir == 0
        theta_new = theta;
        omega_new = omega;
        return
    end
    
    % --- Compute stopping distance ---
    % d_stop = w^2 / (2 * a_max)
    d_stop = (omega^2) / (2 * alpha_max_sat);
    
    % --- Decide accelerate or brake ---
    if abs(err) > d_stop
        % Accelerate toward target
        alpha_sat = dir * alpha_max_sat;
    else
        % Brake
        alpha_sat = -sign(omega) * alpha_max_sat;
    end
    
    % --- Integrate ---
    omega_new = omega + (alpha_sat * dt);
    theta_new = theta + (omega_new * dt);
    
    if (omega_new > omega_max)
        omega_new = omega_max;
    end
end