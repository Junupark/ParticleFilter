function tests = TestParticleFilter
    tests = functiontests(localfunctions);
end

function setupOnce(test_case)
    % for teardown, restoration of original condition
    original_path = path;
    repo_root = fileparts(mfilename('fullpath'));
    addpath(repo_root);

    rng(7);

    dimX = 2;
    dimZ = 1;
    Np = 200;
    A = [1, 0.1; 0, 1];
    H = [1, 0];
    Q = 0.03 * eye(dimX);
    R = 0.1;
    transition_fcn = @(x) A * x;
    measurement_fcn = @(x) H * x;
    likelihood_fcn = @(PF, z, idx) mvnpdf(z, measurement_fcn(PF.particles(:, idx)), R);

    pf = ParticleFilter( ...
        Np, transition_fcn, measurement_fcn, likelihood_fcn, dimX, dimZ, Q, ...
        'x0', [0; 1], ...
        'P0', eye(dimX), ...
        'name', "test_pf");

    test_case.TestData.original_path = original_path;
    test_case.TestData.repo_root = repo_root;
    test_case.TestData.pf = pf;
    test_case.TestData.measurement_fcn = measurement_fcn;
end

function teardownOnce(test_case)
    path(test_case.TestData.original_path);
end

function testConstructAndStateSize(test_case)
    pf = test_case.TestData.pf;
    verifySize(test_case, pf.particles, [pf.dimX, pf.Np]);
    verifySize(test_case, pf.weights, [1, pf.Np]);
    verifyEqual(test_case, sum(pf.weights), 1, 'AbsTol', 1e-12);
end

function testPredictUpdateAndMoments(test_case)
    pf = test_case.TestData.pf.Copy();

    pf.Predict();
    z = 0.25;
    [diverged, varianceW, entropyW, normalization_constant] = pf.UpdateMeasurement(z);
    [mmse, covariance] = pf.MMSE();

    verifyClass(test_case, diverged, 'logical');
    verifyTrue(test_case, isfinite(varianceW));
    verifyTrue(test_case, isfinite(entropyW));
    verifyTrue(test_case, normalization_constant > 0);
    verifySize(test_case, mmse, [pf.dimX, 1]);
    verifySize(test_case, covariance, [pf.dimX, pf.dimX]);
    verifyEqual(test_case, sum(pf.weights), 1, 'AbsTol', 1e-12);
end

function testResampleSystematic(test_case)
    pf = test_case.TestData.pf.Copy();

    pf.Predict();
    pf.UpdateMeasurement(0.1);
    [ess, resampled] = pf.Resample('method', 'systematic', 'ratioNeff', 1);

    verifyLessThanOrEqual(test_case, ess, pf.Np);
    verifyTrue(test_case, resampled);
    verifyEqual(test_case, pf.weights, ones(1, pf.Np) / pf.Np, 'AbsTol', 1e-12);
end

function testProposalSamplerHook(test_case)
    pf = test_case.TestData.pf.Copy();

    pf.Predict('ProposalSampler', @PriorProposal);

    z = test_case.TestData.measurement_fcn([0.2; 1.0]);
    pf.UpdateMeasurement(z);

    verifyEqual(test_case, sum(pf.weights), 1, 'AbsTol', 1e-12);
end
