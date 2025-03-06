function dState = DDE(t, State, Robot, FrictionType, StateDelayed, Controller)
    dq = State(3:4);    
    % FrictionType = Stribeck or Sigmoid
    if FrictionType == "Stribeck"
        friction=@(dq)get_StribeckFriction(dq, Robot);
    elseif FrictionType == "Sigmoid"
        friction=@(dq)get_SigmoidFriction(dq, Robot);
    else
        error("Friction type must be Stribeck or Sigmoid!!!");
    end
    Torque = Controller(StateDelayed);
    ddq = get_MassMatrix(State, Robot) \ (Torque - get_CoriolisVector(State,Robot) - ...
        get_GravityVector(State, Robot) - friction(dq));
    dState = [dq; ddq];
end
