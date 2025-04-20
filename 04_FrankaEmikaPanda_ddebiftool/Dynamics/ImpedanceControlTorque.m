function Torque = ImpedanceControlTorque(Kp, Kd, State, Target)
    q = State(1:7, 1);
    dq = State(8:14, 1);
    alphaMatrix = get_MassMatrix(q);
    betaVector = get_CoriolisVector(q, dq) + get_GravityVector(q) + get_FrictionTorque(dq);
    tauPrime = Kp * (q - Target) + Kd * dq;
    Torque = alphaMatrix * tauPrime + betaVector;
end
