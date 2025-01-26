function dState = ODE(t, State, Torque)
    q = State(1:2, :);
    dq = State(3:4, :);
    ddq = get_MassMatrix(q) \ (Torque - get_CoriolisVector(q, dq) - get_GravityVector(q) - get_FrictionTorque(dq));
    dState = [dq; ddq];
end
