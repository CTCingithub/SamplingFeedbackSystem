function sol = SampledControlSolver(SystemODE, InputFun, SolverOpt)
    % SystemODE is a function with the form SystemODE(t, State, Input)
    % InputFun is a function with the form InputFun(State)
    % SolverOpt is a struct, with 
    %   SolverOpt.initial_condition
    %   SolverOpt.end_time
    %   SolverOpt.sampling_time
    %   SolverOpt.ratio
    %   SolverOpt.ode_setting
    % The return sol is a struct, with
    %   sol.ts
    %   sol.dim
    %   sol.sampling_moments
    %   sol.duration
    %   sol.state_history
    %   sol.state_sampled
    %   sol.input_history
    
    % Sampling Time
    sol.ts = SolverOpt.sampling_time;
    
    % Dimension is the dimension of 2nd order Dynamical Equations
    sol.dim = length(SolverOpt.initial_condition) / 2;
    
    % Time step for recording Time History
    if SolverOpt.ratio < 0.1 && SolverOpt.ratio > 0
        sol.dt = SolverOpt.ratio * SolverOpt.sampling_time;
    else
        sol.dt = 0.1 * SolverOpt.sampling_time;
    end
        
    % Sampling Moments
    sol.sampling_moments = 0:sol.ts:SolverOpt.end_time;
    
    % Time Duration
    sol.duration = 0:sol.dt:SolverOpt.end_time;
    
    % Initialize Time History Matrixes For example, end time is 1.0s,
    % sampling time is 0.4s, dt is 0.2s, sol.state_history records states
    % at 0, 0.2, 0.4, 0.6, 0.8, 1.0s, sol.input_history records inputs at
    % 0, 0.2, 0.4, 0.6, 0.8, 1.0s, sol.state_sampled records states at 0,
    % 0.4, 0.8s, Noting that floor(1.0/0.4)=2, floor(1.0/0.2)=5, Then
    % sol.state_history's and sol.input_history's widths are
    % floor(t_end/dt)+1=ceil(t_end/dt), sol.state_sampled's width is
    % floor(t_end/t_s).
    sol.state_history = zeros(2 * sol.dim, length(sol.duration));
    sol.state_sampled = zeros(2 * sol.dim, length(sol.sampling_moments));
    sol.input_history = zeros(sol.dim, length(sol.duration));
    sol.state_history(:,1) = SolverOpt.initial_condition;
    sol.state_sampled(:,1) = SolverOpt.initial_condition; 
    sol.input_history(:,1) = InputFun(SolverOpt.initial_condition);
    
    tic
    bar = waitbar(0,"Start Solving");
    % In the example above,there're 2 full intervals and 1 remaining
    % interval,which are [0, 0.2, 0.4], [0.4, 0.6, 0.8] and [0.8, 1.0]. The
    % full intervals' length are (ts/dt)+1 (coded as round(sol.ts/sol.dt)
    % to avoid variable type error), the i-th sampling moment's index is
    % i*(ts/dt)+1, the i-th interval is
    % sol.duration((i-1)*(ts/dt)+1:i*(ts/dt)+1)
    
    % For $t \in [0, t_s]$,input is computed with state at $t = 0$, the 1st column of sol.state_sampled.
    % For $t \in [t_s, 2 t_s]$, input is computed with state at $t = t_s$, the 1st column of sol.state_sampled.
    % For $t \in [2 t_s, 3 t_s]$, input is computed with state at $t = 2 t_s$, the 2nd column of sol.state_sampled.
    % ……
    % For $t \in [k t_s, (k + 1) t_s]$, input is computed with state at $ t = (k - 1) t_s$, the kth column of sol.state_sampled.    
    for i = 1:floor(SolverOpt.end_time / sol.ts)
        sol.state_sampled(:,i) = sol.state_history(:,(i-1)*round(sol.ts/sol.dt)+1);
        if i==1
            Input = InputFun(sol.state_sampled(:,1));
        else
            Input = InputFun(sol.state_sampled(:,i-1));
        end
        sol.input_history(:,(i-1)*round(sol.ts/sol.dt)+1:i*round(sol.ts/sol.dt)+1)=repmat(Input,1,round(sol.ts/sol.dt)+1);
        [~,State]=ode45(@(t,State)SystemODE(t,State,Input),...
            sol.duration((i-1)*round(sol.ts/sol.dt)+1:i*round(sol.ts/sol.dt)+1),...
            sol.state_history(:,(i-1)*round(sol.ts/sol.dt)+1),SolverOpt.ode_setting);
        sol.state_history(:,(i-1)*round(sol.ts/sol.dt)+1:i*round(sol.ts/sol.dt)+1)=State(1:end,:)';
        str=[num2str(100*i/floor(SolverOpt.end_time/sol.ts)),'%',' Solved'];
        waitbar(i/floor(SolverOpt.end_time/sol.ts),bar,str);
    end
    waitbar(1, bar, "Solving Final Interval");
    
    % In the example above,there're 2 full intervals and 1 remaining
    % interval,which are [0, 0.2, 0.4], [0.4, 0.6, 0.8] and [0.8, 1.0]. The
    % remaining interval is
    % sol.duration(floor(t_end/ts)*floor(ts/dt)+1:end).
    if floor(SolverOpt.end_time / sol.ts) ~= ceil(SolverOpt.end_time / sol.ts)
        [~,State]=ode45(@(t,State)SystemODE(t,State,Input),...
            sol.duration(floor(SolverOpt.end_time/sol.ts)*round(sol.ts/sol.dt)+1:end),...
            sol.state_history(:,floor(SolverOpt.end_time/sol.ts)*round(sol.ts/sol.dt)+1),SolverOpt.ode_setting);
        sol.input_history(:,floor(SolverOpt.end_time/sol.ts)*round(sol.ts/sol.dt)+1:end)=repmat(Input,1, ...
            length(sol.duration)-floor(SolverOpt.end_time/sol.ts)*round(sol.ts/sol.dt));
        sol.state_history(:,floor(SolverOpt.end_time/sol.ts)*round(sol.ts/sol.dt)+1:end)=State(1:end,:)';
    end
    delete(bar)
    toc
end

