function [perfect_z, noisy_z] = Radar_meas(target_x, opt)
    arguments
        target_x
        opt.origin_x = zeros(3,1)
        opt.model = 'rb'
        opt.R
    end
    dx = target_x(1:2) - opt.origin_x(1:2);
    switch opt.model
        case 'r' % range
            z = sqrt(dx'*dx);
        case 'b' % bearing
            z = wrapToPi(atan2(dx(2), dx(1)) - opt.origin_x(3));
        case 'rb' % both
            z = [sqrt(dx'*dx);
                 wrapToPi(atan2(dx(2), dx(1)) - opt.origin_x(3))];
    end
    perfect_z = z;
    noisy_z = z + mvnrnd(zeros(size(z)), opt.R)';
end
