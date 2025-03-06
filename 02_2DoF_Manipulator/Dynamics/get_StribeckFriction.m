function tau_F = get_StribeckFriction(qdot, Robot)
    % Parameters from Nonlinear Parameter Identification of the Franka
    % Emika PANDA Robot: A Comparative Analysis of Friction Models
    Fc = Robot.Fc;
    Fv = Robot.Fv;
    Fs = Robot.Fs;
    Vs = Robot.Vs;
    Ks = Robot.Ks;

    qdot_1 = qdot(1);
    qdot_2 = qdot(2);

    f1 = sign(qdot_1) * (Fc(1) + ...
        Fs(1) * exp(-abs(qdot_1) / Vs(1) * Ks(1))) + Fv(1) * qdot_1;
    f2 = sign(qdot_2) * (Fc(2) + ...
        Fs(2) * exp(-abs(qdot_2) / Vs(2) * Ks(2))) + Fv(2) * qdot_2;

    tau_F = [f1; f2];
end
