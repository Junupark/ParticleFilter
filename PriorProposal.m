function [next_state, info] = PriorProposal(PF, idx, opt)
    % Prior proposal: follows PF transition model from particle idx.
    arguments
        PF % {mustBeA(PF, 'ParticleFilter')}
        idx (1,1) {mustBeInteger, mustBePositive}
        opt.Q = PF.Q
    end
    
    state = PF.particles(:, idx);
    draw = mvnrnd(zeros(size(state)), opt.Q)';

    next_state = PF.Ffcn(state) + draw;
    info.p_draw = mvnpdf(draw, zeros(size(state)), opt.Q);
    info.idx_sampledfrom = idx;
end
