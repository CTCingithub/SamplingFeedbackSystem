function tau_F = get_SigmoidFriction(qdot, Robot)
    A = Robot.Psi_1;
    alpha = Robot.Psi_2;
    qdotsign = Robot.Psi_3;

    num_of_joints = length(qdot);
    tau_F = zeros(num_of_joints, 1);

    for i = 1:num_of_joints
        tau_F(i) = -A(i) / (1 + exp(-alpha(i) * qdotsign(i))) + ...
            A(i) / (1 + exp(-alpha(i) * (qdot(i) + qdotsign(i))));
    end

end
