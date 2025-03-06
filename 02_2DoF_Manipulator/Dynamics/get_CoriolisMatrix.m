function C_Mat = get_CoriolisMatrix(State, Robot)
    % Import Robot Parameters
    m2 = Robot.m2;
    l1 = Robot.l1;
    lc2 = Robot.lc2;

    % Extract Joint Angles and Velocities
    theta_2 = State(2);
    v_1 = State(3);
    v_2 = State(4);

    % Obtain Coriolis Matrix elements
    C11 = -l1 .* lc2 .* m2 .* v_2 .* sin(theta_2);
    C12 = -l1 .* lc2 .* m2 .* v_1 .* sin(theta_2) - l1 .* lc2 .* m2 .* v_2 .* sin(theta_2);
    C21 = l1 .* lc2 .* m2 .* v_1 .* sin(theta_2);
    C22 = 0;
    C_Mat = [C11, C12; C21, C22];
end
