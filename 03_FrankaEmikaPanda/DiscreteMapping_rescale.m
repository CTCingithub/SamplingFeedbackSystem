function result = DiscreteMapping_rescale(ts, config)
    addpath("PandaDyn\");
    addpath("utils\");
    
    tic
    disp("Start Computing M_d")
    result.M_d=get_MassMatrix(config.target);
    toc
    disp("Done")
    
    tic
    disp("Start Computing K_G")
    q = sym('q_%d',[config.dimension 1]);
    result.K_G=vpa(subs(jacobian( ...
        get_GravityVector(q),q),q, ...
        config.target),config.precision);
    clear q
    result.K_G=double(result.K_G);
    toc
    disp("Done")
    
    tic
    disp("Start Computing K_f")
    result.K_f=get_LinearizedFrictionMatrix();
    toc
    disp("Done")
    
    tic
    disp("Start Computing A")
    result.A_cal=[zeros(config.dimension),eye(config.dimension); ...
        - (ts^2 * result.K_G)/result.M_d,- (ts * result.K_f)/result.M_d];
    result.B_cal=[zeros(config.dimension);ts^2 * eye(config.dimension)];
    result.A=expm(result.A_cal);
    toc
    disp("Done")
    
    if config.visual_Acal_rank
        rank_data=zeros(config.visual_Acal_rank,1);
        for i=1:length(rank_data)+1
            rank_data(i)=rank(mpower(result.A_cal,i-1));
        end
        scatter(0:length(rank_data)-1,rank_data);
        xlabel("Power $i$","Interpreter","latex","FontSize",15)
        ylabel("Rank","Interpreter","latex","FontSize",15)
        title("Rank of $\mathcal{A}^i$","Interpreter","latex","FontSize",25)
        grid on
    end
    
    tic
    disp("Start Computing B")
    syms tau
    % Compute B matrix with
    % diagnolization or jordan blocks or truncation
    if strcmpi(config.MatExpOpt, "diag")
        [V,D] = eig(result.A_cal);
        integral_part=vpaintegral( ...
            exp(D*tau),tau,0,1, ...
            "RelTol",1e-30,"AbsTol", ...
            1e-30);
        integral_part=double(integral_part);
        result.B=result.A*V*integral_part/V*result.B_cal;
    elseif strcmpi(config.MatExpOpt, "Jordan")
        [P,J] = jordan(result.A_cal);
        integral_part=vpaintegral( ...
            MatExpTrunc(J*tau, ...
            config.trunc_rank),tau, ...
            0,1,"RelTol",1e-30, ...
            "AbsTol",1e-30);
        integral_part=double(integral_part);
        result.B=result.A*P*integral_part/P*result.B_cal;
    elseif strcmpi(config.MatExpOpt, "Trunc")
        integral_part=vpaintegral( ...
            MatExpTrunc(result.A_cal*(-tau),config.trunc_rank), ...
            tau,0,1,"RelTol",1e-30,"AbsTol",1e-30);
        integral_part=double(integral_part);
        result.B=result.A*integral_part*result.B_cal;
    else
        error("Invalid method! Use diag, Jordan or Trunc.");
    end
    clear tau
    toc
    disp("Done")
    
    tic
    disp("Start Compiling W")
    result.W=[result.A,result.B; ...
        -config.K_p+result.K_G, ...
        -config.K_d,zeros(config.dimension)];
    toc
    disp("Done")
    
    rmpath("PandaDyn\");
    rmpath("utils\");
end

