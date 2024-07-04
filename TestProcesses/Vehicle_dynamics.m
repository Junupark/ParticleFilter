function x_next = Vehicle_dynamics(x, u, opt)
    arguments
        x (:, 1)
        u (:, 1)
        opt.V
        opt.model %{mustBePositive, mustBeInteger} = 1
        opt.dt {mustBePositive} = 0.1
    end    
    switch opt.model
        % fixed wing % CV % turn rate % vel first
        case 'fw' % [X,Y,V,Psi]
            dp = opt.V*[cos(x(4)); sin(x(4))];
            x_next = x + opt.dt*[dp; 0; u];
            
        % fixed wing % CV % turn rate % turn first
        case 'fw2' % [X,Y,V,Psi]
            next_th = x(4) + opt.dt*u;
            next_p = x(1:2) + opt.V*opt.dt*[cos(next_th); sin(next_th)];
            x_next = [next_p; 0; next_th];
            
        % multirotor % CV % delta turn
        case 'mr' % [X,Y,V,Psi]
            next_th = x(4) + u;
            next_p = x(1:2) + opt.V*opt.dt*[cos(next_th); sin(next_th)];
            x_next = [next_p; 0; next_th];
            
        % multirotor % delta V
        case 'mr2' % [X,Y,Vx,Vy]
            dp = x(3:4)*opt.dt;
            dv = u;
            x_next = x + [dp; dv];
            
        % multirotor % delta acc
        case 'mr3' % [X,Y,Vx,Vy,Ax,Ay]
            dp = x(3:4)*opt.dt;
            dv = x(5:6)*opt.dt;
            da = u;
            x_next = x + [dp;dv;da];
            
        otherwise
            error('no such model');
    end
end