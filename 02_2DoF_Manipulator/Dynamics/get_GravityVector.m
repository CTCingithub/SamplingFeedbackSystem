function G = get_GravityVector(t, State, Robot)
    % Import Robot Parameters
    m1 = Robot.m1;
    m2 = Robot.m2;
    l1 = Robot.l1;
    lc1 = Robot.lc1;
    lc2 = Robot.lc2;
    g = Robot.g;

    % Extract joint angles
    theta_1 = State(1);
    theta_2 = State(2);

    % Obtain elements
    G1 = g .* l1 .* m2 .* cos(theta_1) + ...
        g .* lc1 .* m1 .* cos(theta_1) + ...
        g .* lc2 .* m2 .* cos(theta_1 + theta_2);
    G2 = g .* lc2 .* m2 .* cos(theta_1 + theta_2);
    G = [G1; G2];
end
