function G = get_GravityVector(State, RobotStruct)
    % Import RobotStruct Parameters
    
    m1 = RobotStruct.m1;
    m2 = RobotStruct.m2;
    l1 = RobotStruct.l1;
    lc1 = RobotStruct.lc1;
    lc2 = RobotStruct.lc2;
    g = RobotStruct.g;

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
