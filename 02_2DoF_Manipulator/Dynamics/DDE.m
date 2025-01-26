function dState = DDE(t, State, StateDelayed, Controller)
    q = State(1:2, :);
    dq = State(3:4, :);
    Torque = Controller(StateDelayed);
    ddq = get_MassMatrix(q) \ (Torque - get_CoriolisVector(q, dq) - get_GravityVector(q) - get_FrictionTorque(dq));
    dState = [dq; ddq];
end
