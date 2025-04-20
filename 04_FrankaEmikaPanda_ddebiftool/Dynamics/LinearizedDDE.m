function dState = LinearizedDDE(t, State, StateDelayed, Controller)
    q = State(1:7, :);
    dq = State(8:14, :);
    Torque = Controller(StateDelayed);
    ddq = get_MassMatrix(q) \ (Torque - get_GravityVector(q_Target) - get_LinearizedFrictionTorque(dq));
    dState = [dq; ddq];
end

