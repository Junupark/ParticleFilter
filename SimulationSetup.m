function sim = SimulationSetup(opt)
    arguments
        opt.n_particles {mustBePositive, mustBeInteger} = 500
        opt.r_Kitagawa {mustBePositive} = 0.1;
        opt.r_Vanderpol {mustBePositive} = [1, 1];
        opt.r_Lorenz {mustBePositive} = 0.1;
        opt.r_tracking {mustBePositive} = [0.2, deg2rad(0.5)];
        
        opt.q_Kitagawa {mustBePositive} = 10;
        opt.q_Vanderpol {mustBePositive} = [3, 3];
        opt.q_Lorenz {mustBePositive} = 0.1;
        opt.q_tracking {mustBePositive} = [0.5, 0.5];
        
        opt.initX_Kitagawa
        opt.initX_Vanderpol
        opt.initX_Lorenz
        opt.initX_tracking
        opt.e0_norm {mustBePositive} = 1;
        
        opt.dt {mustBePositive} = 0.1;
        opt.delta_t {mustBePositive} = 0.5;
        opt.tf {mustBePositive} = 50;
        opt.V {mustBePositive} = 10;
        
        opt.ratioNeff {mustBeInRange(opt.ratioNeff, 0, 1, 'exclusive')} = 0.7;
    end
    addpath("./TestProcesses/");
    sim.len = opt.tf / opt.delta_t;
    assert(sim.len == round(sim.len), 'Wrong delta_t & tf');
    
    
    sim.dimX_Kitagawa = 1;
    sim.dimX_Vanderpol = 2;
    sim.dimX_Lorenz = 3;
    sim.dimX_tracking = 4;
    sim.dimZ_Kitagawa = 1;
    sim.dimZ_Vanderpol = 2;
    sim.dimZ_Lorenz = 3;
    sim.dimZ_tracking = 2;
    
    sim.q_Kitagawa = opt.q_Kitagawa;
    sim.Q_Kitagawa = diag(sim.q_Kitagawa)^2;
    sim.r_Kitagawa = opt.r_Kitagawa;
    sim.R_Kitagawa = diag(sim.r_Kitagawa)^2;
    
    sim.q_Vanderpol = opt.q_Vanderpol;
    sim.Q_Vanderpol = diag(sim.q_Vanderpol)^2;
    sim.r_Vanderpol = opt.r_Vanderpol;
    sim.R_Vanderpol = diag(sim.r_Vanderpol)^2;
    
    sim.q_Lorenz = opt.q_Lorenz;
    sim.Q_Lorenz = diag(sim.q_Lorenz)^2;
    sim.r_Lorenz = opt.r_Lorenz;
    sim.R_Lorenz = diag(sim.r_Lorenz)^2;
    
    sim.q_tracking = opt.q_tracking;
    sim.Q_tracking = diag(sim.q_tracking)^2;
    sim.r_tracking = opt.r_tracking;
    sim.R_tracking = diag(sim.r_tracking)^2;
    
    sim.delta_t = opt.delta_t;
    sim.dt = opt.dt;
    sim.tf = opt.tf;
    
    sim.Np = opt.n_particles;
    sim.ratioNeff = opt.ratioNeff;
    
    
    sim.initX_Kitagawa = opt.initX_Kitagawa;
    sim.initX_Vanderpol = opt.initX_Vanderpol;
    sim.initX_Lorenz = opt.initX_Lorenz;
    
    sim.V = opt.V;
    sim.initX_tracking = opt.initX_tracking;
    sim.initX_tracking(3) = sim.V;
    sim.e0_norm = opt.e0_norm;
    
    sim.X0_Kitagawa = sim.initX_Kitagawa + unit(rand(sim.dimX_Kitagawa,1), 'norm', sim.e0_norm);
    sim.X0_Vanderpol = sim.initX_Vanderpol + unit(rand(sim.dimX_Vanderpol,1), 'norm', sim.e0_norm);
    sim.X0_Lorenz = sim.initX_Lorenz + unit(rand(sim.dimX_Lorenz,1), 'norm', sim.e0_norm);
    sim.X0_tracking = sim.initX_tracking + unit(rand(sim.dimX_tracking,1), 'norm', sim.e0_norm);
end

function e = unit(v, opt)
    arguments
        v (:,1)
        opt.norm = 1
    end
    e = opt.norm*v./norm(v);
end