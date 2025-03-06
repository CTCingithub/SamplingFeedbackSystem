function M = get_MassMatrix(State, Robot)
    % Import Robot Parameters
    m1 = Robot.m1;
    m2 = Robot.m2;
    l1 = Robot.l1;
    lc1 = Robot.lc1;
    lc2 = Robot.lc2;
    Jc1 = Robot.Jc1;
    Jc2 = Robot.Jc2;

    % Extract joint angles
    theta_2 = State(2);

    % Compute mass matrix elements
    M11 = Jc1 + Jc2 + l1 .^ 2 .* m2 + ...
        2 * l1 .* lc2 .* m2 .* cos(theta_2) + ...
        lc1 .^ 2 .* m1 + lc2 .^ 2 .* m2;
    M12 = Jc2 + l1 .* lc2 .* m2 .* cos(theta_2) + lc2 .^ 2 .* m2;
    M22 = Jc2 + lc2 .^ 2 .* m2;
    M = [M11, M12; M12, M22];
end
