function dState = LinearizedODE(t, State, Robot, FrictionType, Torque)
    dq = State(3:4);
    % FrictionType = Stribeck or Sigmoid
    if FrictionType == "Stribeck"
        disp("Unfinished, using Sigmoid.")
        friction=get_SigmoidFrictionLinearized(dq, Robot);
    elseif FrictionType == "Sigmoid"
        friction=get_SigmoidFrictionLinearized(dq, Robot);
    else
        error("Friction type must be Stribeck or Sigmoid!!!");
    end
    ddq = get_MassMatrix(State, Robot) \ (Torque - ...
        get_GravityVector(State, Robot) - friction);
    dState = [dq; ddq];
end

function tau_F = get_SigmoidFrictionLinearized(qdot, Robot)
    A = Robot.Psi_1;
    alpha = Robot.Psi_2;
    qdotsign = Robot.Psi_3;

    num_of_joints = length(qdot);
    tau_F = zeros(num_of_joints, 1);

    for i = 1:num_of_joints
        tau_F(i) = ((A(i) * exp(-qdotsign(i) * alpha(i)) * alpha(i)) /...
            (1 + exp(-qdotsign(i) * alpha(i)))^2) * qdot(i);
    end

end
