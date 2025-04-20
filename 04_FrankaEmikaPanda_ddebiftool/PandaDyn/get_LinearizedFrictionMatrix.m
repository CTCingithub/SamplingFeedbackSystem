function LinearizedFrictionMat = get_LinearizedFrictionMatrix()    
    num_of_joints = 7;
    diagonal_element=zeros(num_of_joints,1);
    
    fp = [...
        0.54615,0.87224,0.64068,1.2794,0.83904,0.30301,0.56489;...
        0,0,0,0,0,0,0;...
        0.039533,0.025882,-0.04607,0.036194,0.026226,-0.021047,0.0035526;...
        5.1181,9.0657,10.136,5.5903,8.3469,17.133,10.336];

    A = fp(1,:);            %\phi_{1}
    k = fp(2,:);
    qdotsign = fp(3,:);     %\phi_{3}
    alpha= fp(4,:);         %\phi_{2}
    
    for i=1:num_of_joints
        diagonal_element(i) = ((A(i)*alpha(i)*exp(-alpha(i)*qdotsign(i)))/...
            (1+exp(-alpha(i)*qdotsign(i)))^2);
    end
    
    LinearizedFrictionMat=diag(diagonal_element);
end

