function [perfect_x, noisy_x] = Kitagawa_process(x, k, opt)
    arguments
        x (1,1)
        k (1,1) {mustBePositive}
        opt.Q {mustBePositive} = 10
    end
    
    perfect_x = 0.5*x + 25*x/(1+x^2) + 8*cos(1.2*k);
    noisy_x = perfect_x + mvnrnd(zeros(size(x)), opt.Q)';
end
