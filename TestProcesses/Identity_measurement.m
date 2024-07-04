function [perfect_z, noisy_z] = Identity_measurement(x, opt)
    arguments
        x (:,1)
        opt.R = 1
    end
    
    perfect_z = x;
    if isscalar(opt.R)
        noisy_z = x + mvnrnd(zeros(size(x)), diag(opt.R*ones(size(x))))';
    else
        noisy_z = x + mvnrnd(zeros(size(x)), opt.R)';
    end
end