function Torque = PDTorque(Kp, Kd, State, Target)
    q = State(1:7, 1);
    dq = State(8:14, 1);
    Torque = -Kp * (q - Target) - Kd * dq + get_GravityVector(q);
end
