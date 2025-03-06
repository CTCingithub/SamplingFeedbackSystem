function dState = ODE(t, State, Robot, FrictionType, Torque)
    dq = State(3:4);
    % FrictionType = Stribeck or Sigmoid
    if FrictionType == "Stribeck"
        friction=get_StribeckFriction(dq, Robot);
    elseif FrictionType == "Sigmoid"
        friction=get_SigmoidFriction(dq, Robot);
    else
        error("Friction type must be Stribeck or Sigmoid!!!");
    end
    ddq = get_MassMatrix(State, Robot) \ (Torque - get_CoriolisVector(State, Robot) - ...
        get_GravityVector(State, Robot) - friction);
    dState = [dq; ddq];
end
