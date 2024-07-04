% Maintainer: Junwoo Park <vividlibra@gmail.com>

classdef ParticleFilter < handle
    properties (SetAccess = private, GetAccess = public)
        % particles: row stack of state(column)
        % weights: row stack of weight(single)
        % {i}th particle: particles(:, i)
        % {i}th weight: weights(i)
        particles
        weights
        
        Ffcn
        Hfcn
        linearizedHfcns
        
        Q
        R
        
        Likelihood
        
        Jp
        Jq
        normalization_constant

        mmse
        covariance
        
        name
    end
    properties (Access = private)
        particles_p
        weights_p
        
        hasPropagated
    end
    properties (SetAccess = private, GetAccess = public)
        Np
        dimX
        dimZ
    end
    methods
        function self = ParticleFilter(Np, TransitionFcn, MeasurementFcn, Likelihood, dimX, dimZ, ProcessNoise, opt)
            arguments
                Np (1,1) {mustBeGreaterThan(Np, 10), mustBeInteger}
                TransitionFcn function_handle
                MeasurementFcn function_handle
                Likelihood function_handle
                dimX (1,1) {mustBePositive, mustBeInteger}
                dimZ (1,1) {mustBePositive, mustBeInteger}
                ProcessNoise
                opt.x0 (1,:) double {mustBeNumeric} = zeros(dimX,1)
                opt.P0 double {mustBeNumeric} = zeros(dimX, dimX)
                opt.x0_exhaustive double
                opt.weights (1,:)
                opt.name string
            end
            self.Np = Np;
            self.dimX = dimX;
            self.dimZ = dimZ;
            
            if ~isfield(opt, 'weights')
                self.weights = ones(1,self.Np)/self.Np;
                self.weights_p = ones(1,self.Np)/self.Np;
            else
                assert(length(opt.weights) == self.Np);
                if abs(1 - sum(opt.weights)) > 1e-15
                    warning('Inappropriate weight sum: Assigned after normalized');
                end
                self.weights = opt.weights./sum(opt.weights);
                self.weights_p = self.weights;
            end
            self.normalization_constant = 1;
            
            assert(length(opt.x0)==self.dimX && all(size(opt.P0)==[self.dimX self.dimX]), 'Inappropriate format');
            if isfield(opt, 'x0_exhaustive')
                assert(all(size(opt.x0_exhaustive) == [self.dimX, self.Np]), 'Inappropriate Initial Particles')
                self.particles = opt.x0_exhaustive;
            else
                self.particles = self.DrawNSamplesMVN('mu', opt.x0, 'cov', opt.P0);
            end
            self.particles_p = self.particles;
            self.Jp = zeros(self.Np, 1);
            self.Jq = zeros(self.Np, 1);
            self.Q = ProcessNoise;
            
            assert(all(size(TransitionFcn(self.particles(:, randi(self.Np, 1))))==[self.dimX, 1]), 'Inappropriate Transition Fcn');
            assert(all(size(MeasurementFcn(self.particles(:, randi(self.Np, 1))))==[self.dimZ, 1]), 'Inappropriate Measurement Fcn');
            
            self.Ffcn = TransitionFcn;
            self.Hfcn = MeasurementFcn;
            self.Likelihood = Likelihood;
            self.linearizedHfcns = zeros(self.dimZ, self.dimX, self.Np);
            self.hasPropagated = false;
            
            if isfield(opt, 'name')
                self.name = opt.name;
            else
                self.name = "temporary";
            end
        end
        
        
        
        function self = Predict(self, opt)
            arguments
                self {mustBeA(self, 'ParticleFilter')}
                opt.Q
                opt.TransitionFcn function_handle
                opt.ProposalSampler function_handle
                opt.dry_run logical = false
                opt.verbose logical = false
            end
            if isfield(opt, 'Q')
                assert(all(size(opt.Q) == size(self.Q)))
                self.Q = opt.Q;
            end
            if isfield(opt, 'TransitionFcn')
                self.Ffcn = opt.TransitionFcn;
            end
            
            % Keep up with the latest, Mostly not the case
            if self.hasPropagated; KeepUp(self); end
            
            % <Bootstrap PF>
            if ~isfield(opt, 'ProposalSampler') 
                v = self.DrawNSamplesMVN('cov', self.Q);
                for i=1:self.Np
                    self.particles_p(:,i) = self.Ffcn(self.particles(:,i)) + v(:,i);
                end
                self.normalization_constant = 1;

            % <Proposal PF>
            else
                weights_temp = zeros(size(self.weights));
                for i=1:self.Np
                    [self.particles_p(:,i), info] = opt.ProposalSampler(self, i);
                    e = self.particles_p(:,i) - self.Ffcn(self.particles(:,info.idx_sampledfrom));
                    self.Jp(i) = mvnpdf(e, zeros(self.dimX,1), self.Q);
                    self.Jq(i) = info.p_draw;
                    weights_temp(i) = self.weights(info.idx_sampledfrom)*self.Jp(i)/self.Jq(i);
                    
                    if isfield(info, 'H')
                        self.linearizedHfcns(:,:,i) = info.H;
                    end
                end
                % turn specific NANs into 0
                weights_temp(isnan(weights_temp)) = 0;
                
                % if all was NAN(=0)
                if ~any(weights_temp)
                    self.weights = ones(1,self.Np)/self.Np;
                    self.normalization_constant = 1;
                    if opt.verbose
                        warning('! Filter<%s> RESET @ Proposal', self.name);
                    end
                else
                    self.normalization_constant = self.NormalizeWeights('from', weights_temp);
                end
            end

            if ~opt.dry_run
                self.hasPropagated = true;
            end
        end
        
        function [diverged, varianceW, entropyW, normalization_constant] = UpdateMeasurement(self, z, opt)
            arguments
                self {mustBeA(self, 'ParticleFilter')}
                z
                opt.Likelihood function_handle = self.Likelihood
                opt.verbose = false
                opt.divergence_threshold = 1e-20
            end
            
            % time-propagated case(mostly the case)
            if self.hasPropagated
                self.KeepUp();
                self.hasPropagated = false;
            end
            
            for i=1:self.Np
                self.weights(i) = self.weights(i) * opt.Likelihood(self, z, i);
            end
            
            diverged = sum(self.weights) < opt.divergence_threshold;
            if diverged && opt.verbose
                warning('! Filter<%s> DIVERGED.\tsum(w)<1e%d', self.name, log10(opt.divergence_threshold));
            end
            
            normalization_constant = self.normalization_constant*self.NormalizeWeights('dry_run',true);
            self.normalization_constant = normalization_constant;
            
            if ~any(self.weights) % if all weight was ~0(negligible): Reset
                self.weights = ones(1,self.Np)/self.Np;
                if opt.verbose
                    warning('! Filter<%s> RESET @ Update', self.name);
                end
            else
                self.NormalizeWeights();
            end
            
            varianceW = var(self.weights);
            entropyW = Entropy(self.weights, true);
        end
        
        
        
        function KeepUp(self)
            self.particles(:) = self.particles_p(:);
        end
        function CopyParticlesFromIdx(self, idx)
            if self.hasPropagated
                self.particles = self.particles_p(:,idx);
                self.hasPropagated = false;
            else
                self.particles = self.particles(:,idx);
            end
        end
        function CopyParticlesFrom(self, particles)
            assert(all(size(particles) == [self.dimX, self.Np]))
            self.particles = particles;
        end
        
        function normalization_constant = NormalizeWeights(self, opt)
            arguments
                self
                opt.from (1,:)
                opt.dry_run logical = false
            end
            if isfield(opt, 'from')
                assert(length(opt.from)==self.Np)
                normalization_constant = sum(opt.from);
                if ~opt.dry_run
                    self.weights = opt.from./normalization_constant;
                end
            else
                normalization_constant = sum(self.weights);
                if ~opt.dry_run
                    self.weights = self.weights./normalization_constant;
                end
            end
        end
        
        function [ESS, resampled] = Resample(self, opt)
            arguments
                self {mustBeA(self, 'ParticleFilter')}
                opt.method (1,:) char ... 
                    {mustBeMember(opt.method, {'systematic','stratified','multinomial','residual','autoregulation'})} ...
                    = 'systematic'
                opt.ratioNeff double ...
                    {mustBeInRange(opt.ratioNeff, 0, 1, 'inclusive')} = 1
                opt.rougheningQ double = zeros(self.dimX)
                opt.dr_charging_criteria char {mustBeMember(opt.dr_charging_criteria, {'mean','median','sum_half'})} = 'mean'
                opt.dr_charging_formula char {mustBeMember(opt.dr_charging_formula, {'prop','inv_m_log', 'exp', 'inv', 'm_log'})} = 'prop'
                opt.dr_kernel (1,:) char {mustBeMember(opt.dr_kernel, {'uniform', 'triangular', 'epanechnikov', 'normal'})} = 'epanechnikov'
                opt.dr_bw = 0;
                opt.dr_k function_handle = @Gate
                opt.verbose = false
            end
            
            % Check effective particles
            ESS = self.EffectiveParticles();
            if ESS > self.Np * opt.ratioNeff
                if opt.verbose
                    warning("! Filter<%s> NOT-RESAMPLED", self.name);
                end
                resampled = false;
                return;
            end
            
            % Pre-resample processes for particular PFs
            regularized = contains(self.name, 'regularized', 'IgnoreCase', true);
            if regularized
                self.MMSE(); % to ensure we have updated (empirical) covariance
                D = chol(self.covariance + 1e-17*eye(2), 'lower'); % to ensure the positive definite (not semi-) ness
            end

            % Resample accordingly
            switch opt.method
                case 'systematic'
                    [idx, self.weights] = Resample.systematic(self.weights);
                    self.CopyParticlesFromIdx(idx);
                case 'stratified'
                    [idx, self.weights] = Resample.stratified(self.weights);
                    self.CopyParticlesFromIdx(idx);
                case 'multinomial'
                    [idx, self.weights] = Resample.multinomial(self.weights);
                    self.CopyParticlesFromIdx(idx);
                case 'residual'
                    [idx, self.weights] = Resample.residual(self.weights);
                    self.CopyParticlesFromIdx(idx);
                case 'autoregulation'
                    [particles_, self.weights] = ...
                        DeterministicResample.Synchronization(self, ...
                            'charging_criteria', opt.dr_charging_criteria, 'charging_formula', opt.dr_charging_formula, ...
                            'kernel', opt.dr_kernel, 'k', opt.dr_k, 'bw', opt.dr_bw);
                    self.CopyParticlesFrom(particles_);
                otherwise
                    error('not supported');
            end

            % Post-resample processes for particular PFs
            if regularized % using Gaussian kernel
                p = 1/self.dimX+4;
                A = (4/(self.dimX+2))^p;
                h = A*(self.Np^(-p));
                self.particles = self.particles + h*D*randn(self.dimX,self.Np);
            end
            if any(opt.rougheningQ(:) > 0)
                self.particles = self.particles + self.DrawNSamplesMVN('cov', opt.rougheningQ);
            end
            resampled = true;
        end
        
        function Neff = EffectiveParticles(self)
            Neff = 1/sum(self.weights.^2);
        end
        
        function E = Expectation(self, opt)
            arguments
                self
                opt.of function_handle
                opt.lookahead logical = false
            end
            
            if ~isfield(opt, 'of')
                E = self.MMSE(); return;
            else
                fcn = opt.of;
            end
            if opt.lookahead
                particles_ = self.particles_p;
            else
                particles_ = self.particles;
            end
            
            fcn_ret = fcn(particles_(:,1));
            % dim(ret) * dim(interval)
            
            fcn_ret(:, :, self.Np) = 0;
            for i=2:self.Np
                fcn_ret(:,:,i) = fcn(particles_(:,i));
            end
            E = sum(fcn_ret .* reshape(self.weights, 1, 1, self.Np), 3);
        end
        function [mmse, covariance] = MMSE(self, opt)
            arguments
                self
                opt.lookahead logical = false
            end
            
            if opt.lookahead || self.hasPropagated
                mmse = self.particles_p * self.weights';
                deviation = self.particles_p - mmse;
            else
                mmse = self.particles * self.weights';
                deviation = self.particles - mmse;
            end
            covariance = deviation.*self.weights*deviation';
            self.mmse = mmse;
            self.covariance = covariance;
        end
        
        function samples = DrawNSamplesMVN(self, opt)
            arguments
                self {mustBeA(self, 'ParticleFilter')}
                opt.mu (1,:) = zeros(self.dimX,1)
                opt.cov = zeros(self.dimX)
                opt.N = self.Np
                opt.overwrite logical = false
            end
            assert(all(size(opt.cov)==[length(opt.mu), length(opt.mu)]), 'dimension mismatch');
            samples = mvnrnd(opt.mu, opt.cov, opt.N)';
            if opt.overwrite
                self.CopyParticlesFrom(samples);
                self.weights = ones(1,opt.N)/opt.N;
            end
        end
        
        function Plot(self, opt, opt_plot)
            arguments
                self {mustBeA(self, 'ParticleFilter')}
                opt.propagated logical = false
                opt_plot.?matlab.graphics.chart.primitive.Line
            end
            opt_plot = namedargs2cell(opt_plot);
            if opt.propagated
                particle_plot = self.particles_p;
            else
                particle_plot = self.particles;
            end
            h = scatter(particle_plot(2,:), particle_plot(1,:), 250+400*self.weights, opt_plot{:});
            grid(h.Parent, 'on');
        end
        
        function [P, W, N] = GetEnsemble(self, opt)
            arguments
                self
                opt.lookahead logical = false
            end
            
            if opt.lookahead || self.hasPropagated
                P = self.particles_p;
            else
                P = self.particles;
            end
            
            W = self.weights;
            N = self.Np;
        end
        
        function new_object = Copy(self, opt)
            arguments
                self
                opt.lookahead logical = false
            end
            [P, W, N] = GetEnsemble(self, 'lookahead', opt.lookahead);
            new_object = ParticleFilter(N, self.Ffcn, self.Hfcn, self.Likelihood, self.dimX, self.dimZ, self.Q, ...
                                            'x0_exhaustive', P, 'weights', W);
        end
    end
    
    methods (Static = true)
        function E = Expectation_(particles, weights, opt)
            arguments
                particles
                weights
                opt.of function_handle
            end
%             assert(size(particles, 2) == length(weights));
            if ~isfield(opt, 'of')
                E = ParticleFilter.MMSE_(particles, weights); return;
            else
                fcn = opt.of;
            end
            N = length(weights);
            
            fcn_ret = fcn(particles(:,1));
            % dim(ret) * dim(interval)
            
            fcn_ret(:, :, N) = 0;
            for i=2:N
                fcn_ret(:,:,i) = fcn(particles(:,i));
            end
            E = sum(fcn_ret .* reshape(weights, 1, 1, N), 3);
        end
        function [mmse, covariance] = MMSE_(particles, weights)
            assert(size(particles, 2) == length(weights));
            mmse = particles * weights';
            deviation = particles - mmse;
            covariance = deviation.*weights*deviation';
        end
    end
end 