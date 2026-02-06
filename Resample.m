% Maintainer: Junwoo Park <vividlibra@gmail.com>

% R1: https://bisite.usal.es/archivos/resampling_methods_for_particle_filtering_classification_implementation_and_strategies.pdf
% R2: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4378824
% R3: https://uk.mathworks.com/matlabcentral/fileexchange/24968-resampling-methods-for-particle-filtering
classdef Resample
    methods (Static = true)
        function [idx, weights] = systematic(weights)
            [N, CDF] = PreProcess(weights);
            if all(isnan(CDF))
                idx = randi(N, 1, N);
            else
                U = ((0:N-1)+rand)./N;
                idx = Resample.InverseTransformSampling(CDF, U);
            end
            weights = EqualWeights(N);
        end
        
        function [idx, weights] = stratified(weights)
            [N, CDF] = PreProcess(weights);
            if all(isnan(CDF))
                idx = randi(N, 1, N);
            else
                U = ((0:N-1)+rand(1,N))./N;
                idx = Resample.InverseTransformSampling(CDF, U);
            end
            weights = EqualWeights(N);
        end
        
        function [idx, weights] = multinomial(weights)
            [N, CDF] = PreProcess(weights);
            if all(isnan(CDF))
                idx = randi(N, 1, N);
            else
                U = sort(rand(1,N));
                idx = Resample.InverseTransformSampling(CDF, U);
            end
            weights = EqualWeights(N);
        end
        
        function [idx, weights] = residual(weights)
            N = length(weights);
            if all(isnan(weights))
                idx = randi(N, 1, N);
            else
                weights_scale = weights.*N;
                weights_floor = floor(weights_scale);
                Nr = N - sum(weights_floor);
                weights_residual = (weights_scale - weights_floor)/Nr;

                % Deterministic part
                idx = zeros(1,N); i=1;
                for j=1:N
                    for k=1:weights_floor(j)
                        idx(i)=j;
                        i=i+1;
                    end
                end

                % Stochastic part
                U = sort(rand(1,Nr)); % multinomial sampling
                [~, Cdf] = PreProcess(weights_residual);
                idx(i:end) = Resample.InverseTransformSampling(Cdf, U);
            end
            weights = EqualWeights(N);
        end
    end
    methods (Static = true, Access = private)
        function idx = InverseTransformSampling(Cdf, U)
            N = length(U);
            idx = zeros(1,N);
            j = find(Cdf>0, 1); % filter-out ZERO weight 
            for i = 1:N
                while Cdf(j) < U(i)
                    j = j + 1;
                end
                idx(i) = j;
            end
        end
    end
end

function [N, CDF] = PreProcess(weights)
    N = length(weights);
    CDF = cumsum(weights);
end

function weights = EqualWeights(N)
    arguments
        N {mustBeInteger, mustBePositive}
    end
    weights = ones(1,N)/N;
end
