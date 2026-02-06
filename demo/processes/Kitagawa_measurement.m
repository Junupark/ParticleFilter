function [perfect_z, noisy_z] = Kitagawa_measurement(x, opt)
    arguments
        x (1,1)
        opt.R {mustBePositive} = 1
    end
    
    perfect_z = x^2/20;
    noisy_z = perfect_z + mvnrnd(zeros(size(x)), opt.R)';
end
