function [perfect_x, noisy_x] = VanderPol_process(x, k, opt)
    arguments
        x (2,1)
        k (1,1) {mustBePositive}
        opt.dt {mustBePositive} = 0.02
        opt.delta_t {mustBePositive} = 0.2
        opt.mu = 3.5
        opt.Q = eye(2)
    end
    n_prop = opt.delta_t / opt.dt;
    assert(round(n_prop) == n_prop)
    
    f = @(t, x)f_Vanderpol(x, t, opt);
    for i = 1:n_prop
        x = RK4(f, k, x, opt.dt);
%         x = Euler(f, k, x, opt.dt);
    end
        
    perfect_x = x;
    noisy_x = perfect_x + mvnrnd(zeros(size(x)), opt.Q)';
end

function f = f_Vanderpol(x, t, opt)
    dx1dt = x(2);
    dx2dt = opt.mu*(1-x(1)^2)*x(2) - x(1);
    
    f = [dx1dt; dx2dt];
end