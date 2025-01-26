function C_Vec = get_CoriolisVector(t, State, Robot)
    % Import Robot Parameters
    m2 = Robot.m2;
    l1 = Robot.l1;
    lc2 = Robot.lc2;

    % Extract Joint Angles and Velocities
    theta_2 = State(2);
    v_1 = State(3);
    v_2 = State(4);

    % Obtain Coriolis Vector elements
    C1 = l1 .* lc2 .* m2 .* v_2 .* (-2 * v_1 - v_2) .* sin(theta_2);
    C2 = l1 .* lc2 .* m2 .* v_1 .^ 2 .* sin(theta_2);
    C_Vec = [C1; C2];
end
