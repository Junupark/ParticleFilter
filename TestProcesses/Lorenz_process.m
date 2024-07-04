function [perfect_x, noisy_x] = Lorenz_process(x, k, opt)
    arguments
        x (3,1)
        k (1,1) {mustBePositive}
        opt.dt {mustBePositive} = 0.005
        opt.delta_t {mustBePositive} = 0.02
        opt.sigma = 10
        opt.beta = 8/3
        opt.rho = 28
        opt.Q {mustBePositive} = 2
    end
    n_prop = opt.delta_t / opt.dt;
    assert(round(n_prop) == n_prop)
    
    f = @(t, x)f_Lorenz(x, t, opt);
    for i = 1:n_prop
        x = RK4(f, k, x, opt.dt);
%         x = Euler(f, k, x, opt.dt);
    end
    
    perfect_x = x;
    noisy_x = perfect_x + mvnrnd(zeros(size(perfect_x)), opt.Q)';
end

function f = f_Lorenz(x, t, opt)
    dxdt = opt.sigma*(x(2) - x(1));
    dydt = x(1)*(opt.rho - x(3)) - x(2);
    dzdt = x(1)*x(2) - opt.beta*x(3);
    
    f = [dxdt; dydt; dzdt];
end