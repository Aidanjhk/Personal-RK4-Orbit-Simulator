function omega_dot_new = satAttitudeControl(angularVelSat,angularVelRW,angularAccRW, satInertia,rwInertia, referenceDirections, commandTorque, externalTorque)
 
    W = rwInertia;
    omega_dot_new = -(satInertia^-1)*cross(angularVelSat,satInertia*angularVelSat + W*angularVelRW) + (satInertia^-1)*(externalTorque'-W*angularAccRW')
end