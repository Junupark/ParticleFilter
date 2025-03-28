function [next_state, info] = PriorProposal(state, f, z, opt)
    % f(x_[k-1])
    arguments
        state
        f function_handle
        z 
        opt.Q
    end
    
    draw = mvnrnd(zeros(size(state)), opt.Q);
    
    next_state = f(state) + draw;
    info.p_draw = mvnpdf(draw, zeros(size(state)), opt.Q);
end
