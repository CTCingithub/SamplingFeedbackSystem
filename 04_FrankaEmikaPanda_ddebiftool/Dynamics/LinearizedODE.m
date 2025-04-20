function dState = LinearizedODE(t, State, q_Target, Torque)
    q = State(1:7, :);
    dq = State(8:14, :);
    ddq = get_MassMatrix(q_Target) \ (Torque - get_GravityVector(q_Target) - get_LinearizedFrictionTorque(dq));
    dState = [dq; ddq];
end
