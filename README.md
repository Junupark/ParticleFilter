# ParticleFilter

A MATLAB implementation of modular Particle Filter

## Usage
Sample: [main.m](./main.m)

### 1. Particle Filter class
- input:
  - Number of particles (N)
  - Prediction function
  - Measurement function
  - Likelihood function
  - dimension of state vector (dim(x))
  - dimension of measuremente vector (dim(z))
  - prediction noise covariance (Gaussian), you can use your own sampler (possibly non-Gaussian) at prediction time
  - Keyword Arguments
    - x0: initial estimate
    - P0: initial covariance (this does not mean that the posterior is Gaussian)
    - x0_exhaustive: initial population of particles (dim(x) times N)
    - weights: initial weight (optional, default: all 1/N)
    - name: string name of the instance (optional)
    - either x0/P0 or x0_exhaustive must be given

### 2. Prediction function

any function: `@(x_{k}) -> x_{k+1}`
- input: a state vector of size dim(x)
- output: next state of size dim(x)

### 3. Measurement function
any function: `@(x_{k}) -> z_{k}`
- input: a state vector of size dim(x)
- output: a measurement vector of size dim(z)

### 4. Likelihood function
any function: `@(PF, z, index) -> p(z_{k} | x_{k}^{index})`
- input: PF instance, measurement, particle index
- output: likelihood value


## Examples 
### Example: Prediction, Measurement, and Likelihood functions
```MATLAB
R = eye(2);
z = rand(2,1);
Prediction = @(x)x;
Measurement = @(x)x; # dim(z) = dim(x)
Likelihood = @(PF, z, index)mvnpdf(z - Measurement(PF.particles(:,index)), R);
% likelihood need not to be Gaussian, it is just an example
```

### Example: Instantiation
```MATLAB
N = 1000
null_process = @(x)x;
Likelihood_Kitagawa = @(PF, z, index)mvnpdf(z, Kitagawa_measurement(PF.particles(:,index)), 1);
Measurement_Kitagawa = @(x)Kitagawa_measurement(x);
pf = ParticleFilter(N, null_process, Measurement_Kitagawa, Likelihood_Kitagawa, 1, 1, 1, 'x0_exhaustive', mvnrnd(1, 5, N)');
```

### Example: Prediction
```MATLAB
% option1: prediction using initial prediction function & prediction covariance
pf.Predict();

% option2: pass custom prediction function
i_sim = 10; % a sample environment variable managed by user, here indicating 10-th time step k
Ffcn_Kitagawa = @(x)Kitagawa_process(x, i_sim);
pf.Predict('TransitionFcn', Ffcn_Kitagawa);

% option3: pass general sampler = importance sampling & weighting
z = rand();
Ffcn_Kitagawa = @(x)Kitagawa_process(x, i_sim); % baseline
sampler = @(PF, idx)PriorProposal(PF.particles(:,idx), Ffcn_Kitagawa, z, 'Q', 1);
pf.Predict('TransitionFcn', Ffcn_Kitagawa, 'ProposalSampler', sampler);
% this example (prior proposal) is practically & statistically the same as option2, as the name suggests
```
Here, `sampler` is any function: `@(PF, index) -> [x_{k+1}^{index}, information about this sample]`
- input: Particle Filter instance, particle index (index-th particle is sampled)
- output: [next state (of size dim(x)), info]
  - info: a structure composed of
    - p_draw: probability that `sampler` sample that next state, i.e., `q(x_{k+1}^{index}|x_{k}^{index})`
    - (many more you want to relay)
- the weights are adjusted according to the importance sampling, i.e., multiplied by `p/q`

You can use this functionality for custom non-standard sampler, such as
- regularized PF (RPF)
- auxiliary PF (APF)
- or any sampler that leverages current measurement `z_{k}`, i.e., the optimal proposal


### Example: Measurement Update
```MATLAB
z = rand();

% option1: weight update using initial likelihood function
pf.UpdateMeasurement(z);

% option2: pass custom likelihood
Likelihood_custom = @(PF, z, idx)mvnpdf(z - PF.particles(:,idx), 5); % different R
pf.UpdateMeasurement(z, 'Likelihood', Likelihood_custom);
```

## Useful methods
- Resampling
```MATLAB
pf.Resample('method', 'systematic' || 'stratified' || 'multinomial' || 'residual') % resampling
pf.Resample('ratioNeff', 0.7) % resampled only when effective sample size (ESS) is less than 0.7N
% if 'ratioNeff` is not given, this method always resample
```
- Statistics (MMSE, Expectation, or Ensemble)
```MATLAB
pf.Expectation('of', arbitrary_function) % returns the expectation of given function w.r.t the current distribution
pf.MMSE() % returns mean & covariance
pf.GetEnsemble() % returns particle, weight, number of particle
```
