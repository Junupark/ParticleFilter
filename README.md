# ParticleFilter

A MATLAB implementation of modular Particle Filter

## Quick Start

See [main.m](./main.m) for a complete example.

```MATLAB
% number of particles
N = 1000;
% Create particle filter with Gaussian initial distribution
pf = ParticleFilter(N, prediction_fcn, measurement_fcn, likelihood_fcn, ...
                   state_dim, measurement_dim, process_noise_cov, ...
                   'x0', initial_state, 'P0', initial_covariance);

% Alternative: Create with custom initial particle distribution
initial_particles = rand(state_dim, N);  % or any custom distribution
pf = ParticleFilter(N, prediction_fcn, measurement_fcn, likelihood_fcn, ...
                   state_dim, measurement_dim, process_noise_cov, ...
                   'x0_exhaustive', initial_particles);

% Prediction step
pf.Predict();

% Measurement update
pf.UpdateMeasurement(measurement);

% Resample adaptively
pf.Resample('ratioNeff', 0.7);
```

## Related Work

**This implementation has been used in the following research. If this package contributes to your work, you may find these papers relevant:**

**Journal Articles:**
```bibtex
@article{park2025contours,
  title={Contours-seeking Proposal Density Particle Filter and Resilient Terrain-referenced Navigation},
  author={Park, Junwoo and Bang, Hyochoong},
  journal={IEEE Transactions on Aerospace and Electronic Systems},
  year={2025},
  doi={10.1109/TAES.2025.3588823},
  note={in press}
}
```

**Conference Papers:**
```bibtex
@inproceedings{park2022balanced,
  title={Balanced cooperative target search of mobile sensor fleet under localization uncertainty},
  author={Park, Junwoo and Lee, Ho-Hyeong and Bang, Hyochoong},
  booktitle={2022 IEEE 61st Conference on Decision and Control (CDC)},
  pages={6626--6631},
  year={2022},
  organization={IEEE},
  address={Cancun, Mexico},
  doi={10.1109/CDC51059.2022.9992977}
}

@inproceedings{park2022evenly,
  title={Evenly weighted particle filter for terrain-referenced navigation using gaussian mixture proposal distribution},
  author={Park, Junwoo and Bang, Hyochoong},
  booktitle={2022 International Conference on Unmanned Aircraft Systems (ICUAS)},
  pages={177--183},
  year={2022},
  organization={IEEE},
  address={Dubrovnik, Croatia},
  doi={10.1109/ICUAS54217.2022.9836197}
}
```

## Core Components

### ParticleFilter Class
**Required inputs:**
- `Np`: Number of particles (integer > 10)
- `TransitionFcn`: State transition function `@(x_{k}) -> x_{k+1}`
- `MeasurementFcn`: Measurement function `@(x_{k}) -> z_{k}`
- `Likelihood`: Likelihood function `@(PF, z, index) -> p(z_{k} | x_{k}^{index})`
- `dimX`: Dimension of state vector (positive integer)
- `dimZ`: Dimension of measurement vector (positive integer)
- `ProcessNoise`: Process noise covariance (Gaussian)
- `opt`: Options struct containing initialization parameters

**Initialization options (in opt struct):**
- `x0`, `P0`: Initial estimate and covariance (required unless x0_exhaustive is provided)
- `x0_exhaustive`: Initial particle population (dimX × Np)
- `weights`: Initial particle weights (optional, default: uniform)
- `name`: String name for the filter instance (optional)

### Key Methods

#### Prediction:
```MATLAB
pf.Predict();  % Using default transition function
pf.Predict('TransitionFcn', custom_fcn);  % Custom transition
pf.Predict('ProposalSampler', sampler);  % Importance sampling with custom proposal
```

**Custom Proposal Samplers:**
The `ProposalSampler` enables advanced particle filtering techniques by allowing custom sampling strategies:

```MATLAB
% Example: Prior proposal (equivalent to standard prediction)
sampler = @(PF, idx) PriorProposal(PF.particles(:,idx), transition_fcn, measurement, 'Q', process_noise);
pf.Predict('ProposalSampler', sampler);

% Custom sampler function signature: @(PF, index) -> [next_state, info]
% - PF: ParticleFilter instance
% - index: particle index being sampled
% - next_state: predicted state (size: state_dim × 1)
% - info.p_draw: proposal probability q(x_{k+1} | x_{k}, z_{k})
```

**Applications:**
- **Optimal proposals** leveraging current measurement `z_{k}`
- **Auxiliary Particle Filter (APF)**, or **Regularized Particle Filter (RPF)**
- **Custom non-Gaussian samplers** for specialized dynamics

#### Measurement Update:
```MATLAB
pf.UpdateMeasurement(z);  % Using default likelihood
pf.UpdateMeasurement(z, 'Likelihood', custom_likelihood);
```

**Custom Likelihood Functions:**
The likelihood function is completely flexible and does not need to be Gaussian:

```MATLAB
% Likelihood function signature: @(PF, z, index) -> p(z | x^{index})
% - PF: ParticleFilter instance containing particles
% - z: current measurement vector
% - index: particle index to evaluate likelihood for

% Gaussian likelihood (common case) - using custom measurement function
l_gaussian = @(PF, z, idx) mvnpdf(z - custom_meas_fcn(PF.particles(:,idx)), R);

% Non-Gaussian examples:
% Laplace distribution for robust estimation
l_laplace = @(PF, z, idx) prod(exp(-abs(z - PF.Hfcn(PF.particles(:,idx)))/sigma));

% Gaussian mixture for multimodality - using custom measurement function
l_gmm = @(PF, z, idx) w1*mvnpdf(z, custom_meas_fcn1(PF.particles(:,idx)), R1) + ...
                      w2*mvnpdf(z, custom_meas_fcn2(PF.particles(:,idx)), R2);
```

**Applications:**
- **Robust estimation** with heavy-tailed distributions (Laplace, Student-t)
- **Terrain-referenced navigation** with custom terrain models
- **Multi-modal measurements** using mixture distributions
- **Outlier rejection** with custom robust likelihoods

#### Resampling:
**Supported Methods**: systematic, stratified, multinomial, and residual
```MATLAB
pf.Resample('method', 'systematic');
pf.Resample('ratioNeff', 0.7);  % Adaptive resampling
```

#### Statistics:
```MATLAB
[mean, cov] = pf.MMSE();  % Mean and covariance
expectation = pf.Expectation('of', @(x) x.^2);  % Custom expectation
[particles, weights, N] = pf.GetEnsemble();  % Raw particles
```
