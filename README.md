# ParticleFilter
tldr; see main.m

## prediction function:
any function satisfying:
- input: a state vector
- output: next state

## likelihood function:
any function satisfying:
- input: PF instance, measurement, particle index
- output: likelihood value

#### Example
```MATLAB
R = eye(2);
z = rand(2,1);
Prediction = @(x)x;
Measurement = @(x)x;
Likelihood = @(PF, z, idx)mvnpdf(z - Measurement(PF.particles(:,idx)), R);
% likelihood need not to be Gaussian, it is just an example
```

## Particle Filter
- input:
  - Number of particles
  - Prediction function
  - Measurement function
  - Likelihood function
  - dimension of state vector
  - dimension of measuremente vector
  - prediction noise covariance (Gaussian), you can use your own sampler (possibly non-Gaussian) at prediction time
  - Keyword Arguments
    - x0: initial estimate
    - P0: initial covariance (this does not mean that the posterior is Gaussian)
    - x0_exhaustive: dx times N initial population of particles 
    - weights: initial weight (optional)
    - name: string name of the instance (optional)
    - either x0/P0 or x0_exhaustive must be given

#### Example
```MATLAB
null_process = @(x)x;
Likelihood_Kitagawa = @(PF, z, idx)mvnpdf(z, Kitagawa_measurement(PF.particles(:,idx)), 1);
pf_Kitagawa_standard = ParticleFilter(250, null_process, @(x)Kitagawa_measurement(x), Likelihood_Kitagawa, 1, 1, 1, 'x0_exhaustive', mvnrnd(1, 5, 250)');
```

### Prediction
```MATLAB
% option1: prediction using initial prediction function & prediction covariance
pf_Kitagawa_standard.Predict();

% option2: pass custom prediction function
Ffcn_Kitagawa = @(x)Kitagawa_process(x, i_sim);
pf_Kitagawa_standard.Predict('TransitionFcn', Ffcn_Kitagawa);

% option3: pass arbitrarily general sampler
% jwpark: I don't think you are going to use this
sampler = "any function satisfying this API: @(PF, idx) -> [next_state, info], sampline idx-th particle"
pf_Kitagawa_standard.Predict('TransitionFcn', Ffcn_Kitagawa, 'ProposalSampler', sampler);
```

### Measurement Update
```MATLAB
z = rand();

% option1: weight update using initial likelihood function
pf_Kitagawa_standard.UpdateMeasurement(z);

% option2: pass custom likelihood
Likelihood_custom = @(PF, z, idx)mvnpdf(z - PF.particles(:,idx), 4);
pf_Kitagawa_standard.UpdateMeasurement(z, 'Likelihood', Likelihood_custom);
```

### Useful APIs
```MATLAB
pf.Resample('method', 'systematic || 'stratified' || 'multinomial' || 'residual') % resampling
pf.Resample('ratioNeff', 0.7) % resampled only when ESS is less than 0.7N
pf.Expectation('of', arbitrary_function) % expectation of given function w.r.t the distribution
pf.MMSE() % mean & covariance
pf.GetEnsemble() % particle, weight, number of particle
```