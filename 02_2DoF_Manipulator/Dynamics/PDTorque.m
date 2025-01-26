function Torque = PDTorque(Kp, Kd, State, Target)
    q = State(1:2, :);
    dq = State(3:4, :);
    Torque = -Kp * (q - Target) - Kd * dq + get_GravityVector(Target);
end
