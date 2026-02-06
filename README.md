# ParticleFilter

A self-contained MATLAB implementation of particle filter for general nonlinear and non-Gaussian state estimation.

## Quick Start (Standalone)

```matlab
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

% Exploit result
[xhat, Phat] = pf.MMSE();
```

**Run demonstration** (see [demo](./demo/RunParticleFilterDemo.m) for comprehensive examples):

```matlab
run(fullfile(pwd, 'demo', 'RunParticleFilterDemo.m'))
```

**Run tests:**

```matlab
runtests
```

## Citation

If this repository contributes to your research, please cite relevant publications below.

Journal Articles:

```bibtex
@article{Park2025ContoursSeeking,
  title={Contours-Seeking Proposal Density Particle Filter and Resilient Terrain-Referenced Navigation},
  author={Park, Junwoo and Bang, Hyochoong},
  year={2025},
  journal={IEEE Transactions on Aerospace and Electronic Systems},
  volume={61},
  number={6},
  pages={15627--15641},
  doi={10.1109/TAES.2025.3588823}
}

@article{Park2024Sparse,
  title={Sparse Instantiation of Bias Nodes for Factor Graph-based Terrain-referenced Navigation},
  author={Park, Junwoo and Bang, Hyochoong},
  journal={International Journal of Control, Automation and Systems},
  year={2024},
  volume={22},
  number={11},
  pages={3364--3376},
  doi={10.1007/s12555-023-0784-x}
}
```

Conference Papers:

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

## Core Configuration and Methods

`ParticleFilter` constructor supports:

- `Np`: number of particles.
- `TransitionFcn`: state transition function handle.
- `MeasurementFcn`: measurement function handle.
- `Likelihood`: likelihood function handle.
- `dimX`, `dimZ`: state and measurement dimensions.
- `ProcessNoise`: process noise covariance.
- `x0`, `P0`: Gaussian initial state and covariance.
- `x0_exhaustive`: user-supplied initial particle matrix (`dimX x Np`).
- `weights`: optional initial normalized weights.
- `name`: optional filter instance name.

### Predict

```matlab
% default propagation with PF transition model
pf.Predict();

% override transition model for this step
pf.Predict('TransitionFcn', custom_transition_fcn);
```

With proposal sampling:

```matlab
% built-in prior proposal (does not require z)
pf.Predict('ProposalSampler', @PriorProposal);
```

Custom proposal samplers must return:

- `next_state`: proposed state for particle `idx`
- `info.p_draw`: proposal density q(x*{k} | x*{k-1}, z\_{k})
- `info.idx_sampledfrom`: ancestor particle index

### Measurement Update

```matlab
% default likelihood (set in constructor)
[diverged, varianceW, entropyW, c] = pf.UpdateMeasurement(z);

% override likelihood for one update
[diverged, varianceW, entropyW, c] = pf.UpdateMeasurement(z, 'Likelihood', custom_likelihood_fcn);
```

Likelihood signature:

```matlab
custom_likelihood_fcn = @(PF, z, idx) p_z_given_x;
```

### Resample

```matlab
% supported methods: systematic, stratified, multinomial, residual
[ess, did_resample] = pf.Resample('method', 'systematic', 'ratioNeff', 0.7);
```

### Statistics and Ensemble Access

```matlab
% posterior mean/covariance
[xhat, Phat] = pf.MMSE();

% generic expectation E[f(x)]
E = pf.Expectation('of', @(x) x.^2);

% raw particles/weights
[P, W, N] = pf.GetEnsemble();
```

## Submodule Use (Recommended)

This repository is intended to be consumed as a submodule inside a parent application.

Assume parent project layout:

```text
MyApp/
  your_app/
  ParticleFilter/   % this repository as submodule
```

From `MyApp`:

```bash
git submodule add https://github.com/Junupark/ParticleFilter.git
git submodule update --init --recursive
```

In MATLAB (from `MyApp`):

```matlab
addpath(fullfile(pwd, 'ParticleFilter'))
% construct and use ParticleFilter directly from your app entrypoint
```
