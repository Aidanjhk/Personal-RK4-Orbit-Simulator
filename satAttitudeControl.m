function omega_dot_new = satAttitudeControl( ...
    angularVelSat, angularVelRW, angularAccRW, ...
    satInertia, rwInertia, referenceDirections, ...
    commandTorque, externalTorque)

    % Wheel axis matrix (3xN)
    W = referenceDirections;

    % ---- enforce column vectors ----
    angularVelSat = angularVelSat(:);
    angularVelRW  = angularVelRW(:);
    angularAccRW  = angularAccRW(:);
    externalTorque = externalTorque(:);

    % ---- make rwInertia usable ----
    if isvector(rwInertia)
        Jw = diag(rwInertia(:));   % Nx1 -> NxN
    else
        Jw = rwInertia;            % already NxN
    end

    % ---- total angular momentum ----
    H = satInertia * angularVelSat + W * (Jw * angularVelRW);

    % ---- reaction torque from wheels ----
    tau_rw = - W * (Jw * angularAccRW);

    % ---- rigid body dynamics ----
    omega_dot_new = satInertia \ ( ...
        externalTorque + tau_rw ...
        - cross(angularVelSat, H) ...
    );
end