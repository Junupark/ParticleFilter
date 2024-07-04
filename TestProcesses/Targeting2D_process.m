function x_next = Targeting2D_process(x, k, opt)
    arguments
        x (:,1)
        k (1,1) {mustBePositive}
        opt.dt = 0.01
        opt.model %{mustBePositive, mustBeInteger} = 1
    end
    
    switch opt.model
        case 'st' % stationary
            x_next = x;
        case 'Vxy' % [X,Y,Vx,Vy]
            x_next = [x(1) + x(3)*opt.dt;
                    x(2) + x(4)*opt.dt;
                    x(3); 
                    x(4)];
        case 'CV'
            psi = wrapToPi(x(4));
            dp = x(3)*opt.dt*[cos(psi); sin(psi)];
            x_next = x + [dp; 0; 0];
            
        case 'CT'
            psi = wrapToPi(x(4));
            dp = x(3)*opt.dt*[cos(psi); sin(psi)];
            x_next = x + [dp; 0; x(5); 0];
            
        otherwise
            error('No such model');
    end
end