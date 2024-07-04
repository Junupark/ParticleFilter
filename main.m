%%
clear;
global sim
sim = SimulationSetup('n_particles', 500, 'ratioNeff', 0.7, ...
                        'q_Kitagawa', 5, 'r_Kitagawa', 1, ...
                        'q_Vanderpol', 0.1*ones(2,1), 'r_Vanderpol', 0.1*ones(2,1), ...
                        'q_Lorenz', 0.5*ones(3,1), 'r_Lorenz', 0.5*ones(3,1), ...
                        'q_tracking', [0.5;0.5;0.25;deg2rad(5)], 'r_tracking', [0.5; deg2rad(0.5)], ...
                        'initX_Kitagawa', 5, ...
                        'initX_Vanderpol', [2;1], ...
                        'initX_Lorenz', 10*rand(3,1), ...
                        'initX_tracking', [100;50;10;0], ...
                        'e0_norm', 0.1, ...
                        'V', 5, ...
                        'tf', 50);

%% Initialize
null_process = @(x)x;
% =========== Kitagawa ===========
Likelihood_Kitagawa = @(PF, z, idx)mvnpdf(z, Kitagawa_measurement(PF.particles(:,idx)), sim.R_Kitagawa);
pf_Kitagawa_standard = ParticleFilter(sim.Np, null_process, @(x)Kitagawa_measurement(x), Likelihood_Kitagawa, sim.dimX_Kitagawa, ...
                                        sim.dimZ_Kitagawa, sim.Q_Kitagawa, 'x0_exhaustive', mvnrnd(sim.X0_Kitagawa, 5, sim.Np)');
% =========== Vanderpol ===========
Likelihood_Vanderpol = @(PF, z, idx)mvnpdf(z, Identity_measurement(PF.particles(:,idx)), sim.R_Vanderpol);
pf_Vanderpol_stadnard = ParticleFilter(sim.Np, null_process, @(x)Identity_measurement(x), Likelihood_Vanderpol, sim.dimX_Vanderpol, ...
                                        sim.dimZ_Vanderpol, sim.Q_Vanderpol, 'x0_exhaustive', mvnrnd(sim.X0_Vanderpol, 1^2*eye(sim.dimX_Vanderpol), sim.Np)');
% =========== Lorenz ===========
Likelihood_Lorenz = @(PF, z, idx)mvnpdf(z, Identity_measurement(PF.particles(:,idx)), sim.R_Lorenz);
pf_Lorenz_stadnard = ParticleFilter(sim.Np, null_process, @(x)Identity_measurement(x), Likelihood_Lorenz, sim.dimX_Lorenz, ...
                                        sim.dimZ_Lorenz, sim.Q_Lorenz, 'x0_exhaustive', mvnrnd(sim.X0_Lorenz, 2.5^2*eye(sim.dimX_Lorenz), sim.Np)');
% =========== Tracking ===========
Likelihood_tracking = @(PF, z, idx)mvnpdf(z, Radar_meas(PF.particles(:,idx), 'model', 'rb', 'R', sim.R_tracking), sim.R_tracking);
pf_tracking_standard = ParticleFilter(sim.Np, null_process, @(x)Radar_meas(x, 'R', sim.R_tracking), Likelihood_tracking, sim.dimX_tracking, ...
                                        sim.dimZ_tracking, sim.Q_tracking, 'x0_exhaustive', mvnrnd(sim.X0_tracking, 4*sim.Q_tracking, sim.Np)' );

trueTraj_Kitagawa = zeros(sim.dimX_Kitagawa, sim.len);
trueTraj_Vanderpol = zeros(sim.dimX_Vanderpol, sim.len);
trueTraj_Lorenz = zeros(sim.dimX_Lorenz, sim.len);
trueTraj_tracking = zeros(sim.dimX_tracking, sim.len);
trueTraj_Kitagawa(:, 1) = sim.initX_Kitagawa;
trueTraj_Vanderpol(:, 1) = sim.initX_Vanderpol;
trueTraj_Lorenz(:, 1) = sim.initX_Lorenz;
trueTraj_tracking(:, 1) = sim.initX_tracking;

estTraj_Kitagawa = zeros(sim.dimX_Kitagawa, sim.len);
estTraj_Vanderpol = zeros(sim.dimX_Vanderpol, sim.len);
estTraj_Lorenz = zeros(sim.dimX_Lorenz, sim.len);
estTraj_tracking = zeros(sim.dimX_tracking, sim.len);

estTraj_Kitagawa(:, 1) = pf_Kitagawa_standard.MMSE();
estTraj_Vanderpol(:, 1) = pf_Vanderpol_stadnard.MMSE();
estTraj_Lorenz(:, 1) = pf_Lorenz_stadnard.MMSE();
estTraj_tracking(:, 1) = pf_tracking_standard.MMSE();

variance_Kitagawa = zeros(1, sim.len);
variance_Vanderpol = zeros(1, sim.len);
variance_Lorenz = zeros(1, sim.len);
variance_tracking = zeros(1, sim.len);

ess_Kitagawa = zeros(1, sim.len);
ess_Vanderpol = zeros(1, sim.len);
ess_Lorenz = zeros(1, sim.len);
ess_tracking = zeros(1, sim.len);

Z_Kitagawa = zeros(sim.dimZ_Kitagawa, sim.len);
Z_Vanderpol = zeros(sim.dimZ_Vanderpol, sim.len);
Z_Lorenz = zeros(sim.dimZ_Lorenz, sim.len);
Z_tracking = zeros(sim.dimZ_tracking, sim.len);


colors = [0,0,1;1,0,0];
markers = {'d', 's'};
%% ========================Tracking Example========================
disp("Tracking example");
% Simulate True Trajectory & Measurement
for i_sim=2:sim.len
    trueTraj_tracking(:, i_sim) = Vehicle_dynamics(trueTraj_tracking(:, i_sim-1), deg2rad(30*(2*rand-1.2)), 'model', 'fw', 'V', sim.V); % slightly biased
    [~, Z_tracking(:, i_sim)]   = Radar_meas(trueTraj_tracking(1:2, i_sim), 'model', 'rb', 'R', sim.R_tracking);
end

% PF
for i_sim=2:sim.len
    % Propagation
    Ffcn_tracking = @(x)Targeting2D_process(x, i_sim, 'dt', sim.dt, 'model', 'CV');
    pf_tracking_standard.Predict('TransitionFcn', Ffcn_tracking);
    
    % Measurement Update
    [~, variance_tracking(i_sim), ~] = pf_tracking_standard.UpdateMeasurement(Z_tracking(:, i_sim));
    % Resample
    ess_tracking(i_sim) = pf_tracking_standard.Resample('method', 'systematic', 'ratioNeff', sim.ratioNeff);

    estTraj_tracking(:, i_sim) = pf_tracking_standard.MMSE();
    fprintf("%dth step done\n",i_sim);
end
% plot trajectory
fig_tracking = figure('Name', 'Tracking Example');
h_true = plot(trueTraj_tracking(1,1:sim.len),trueTraj_tracking(2,1:sim.len), 'k', 'LineWidth', 3); hold on; axis equal; grid on;
h_meas = plot(Z_tracking(1, 2:sim.len).*cos(Z_tracking(2, 2:sim.len)), Z_tracking(1, 2:sim.len).*sin(Z_tracking(2, 2:sim.len)), 'm:',...
                'LineWidth', 1.5); 
h_pf = plot(estTraj_tracking(1, 1:sim.len),estTraj_tracking(2, 1:sim.len), 'Color', colors(1, :), 'Marker', markers{1}, 'MarkerIndices', 1:10:sim.len); hold on;
legend([h_true, h_meas, h_pf], 'True', 'h^{-1}(y)', 'Standard');
xlabel('x (m)'); ylabel('y (m)');

% Error
figure;
for j = 1:4
    subplot(4,1,j);
    plot([0, sim.len], zeros(1,2), 'k', 'LineWidth', 3); hold on; axis equal; grid on;
    plot(estTraj_tracking(j, 1:sim.len) - trueTraj_tracking(j,1:sim.len), 'Color', colors(1, :), 'LineWidth', 2, 'Marker', markers{1}, 'MarkerIndices', 1:10:sim.len);
end
% Resample Triggered
figure;
subplot(2,1,1);
plot(1:sim.len, ess_tracking);
title('ess');
% Weight Variance
subplot(2,1,2);
plot(1:sim.len, variance_tracking);
title('weight variance');


%% ========================KITAGAWA Example========================
disp("Kitagawa example");
% Simulate True Trajectory & Measurement
for i_sim=2:sim.len
    trueTraj_Kitagawa(:, i_sim)	= Kitagawa_process(trueTraj_Kitagawa(:, i_sim-1), i_sim);
    [~, Z_Kitagawa(:, i_sim)]   = Kitagawa_measurement(trueTraj_Kitagawa(:, i_sim), 'R', sim.R_Kitagawa);
end

% PF
for i_sim=2:sim.len
    % Propagation
    Ffcn_Kitagawa = @(x)Kitagawa_process(x, i_sim);
    pf_Kitagawa_standard.Predict('TransitionFcn', Ffcn_Kitagawa);
    
    % Measurement Update
    pf_Kitagawa_standard.UpdateMeasurement(Z_Kitagawa(:, i_sim));
    
    % Resample
    ess_Kitagawa(i_sim) = pf_Kitagawa_standard.Resample('method', 'systematic', 'ratioNeff', sim.ratioNeff);

    estTraj_Kitagawa(:, i_sim) = pf_Kitagawa_standard.MMSE();
    disp([num2str(i_sim), ' th step done. ']);
end

figure;
plot(trueTraj_Kitagawa, 'Color', colors(1,:), 'LineWidth', 1); hold on; grid on;
plot(estTraj_Kitagawa, 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--');
legend('true', 'estimate');

figure;
plot(ess_Kitagawa);
title('ess');

%% ========================VANDERPOL OSCILLATOR========================
disp("Oscillator example");
% Simulate True Trajectory & Measurement
for i_sim=2:sim.len
    trueTraj_Vanderpol(:, i_sim)    = VanderPol_process(trueTraj_Vanderpol(:, i_sim-1), i_sim);
    [~, Z_Vanderpol(:, i_sim)]      = Identity_measurement(trueTraj_Vanderpol(:, i_sim), 'R', sim.R_Vanderpol);
end

% PF
for i_sim=2:sim.len
    % Propagation
    Ffcn_Vanderpol = @(x)VanderPol_process(x, i_sim);
    pf_Vanderpol_stadnard.Predict('TransitionFcn', Ffcn_Vanderpol);
    
    % Measurement Update
    pf_Vanderpol_stadnard.UpdateMeasurement(Z_Vanderpol(:, i_sim));
    
    % Resample
    pf_Vanderpol_stadnard.Resample('method', 'systematic', 'ratioNeff', sim.ratioNeff);

    estTraj_Vanderpol(:, i_sim) = pf_Vanderpol_stadnard.MMSE();
    disp([num2str(i_sim), ' th step done. ']);
end

figure;
for j=1:2
    subplot(2,1,j);
    plot(trueTraj_Vanderpol(j, :), 'Color', colors(1,:), 'LineWidth', 1); hold on; grid on;
    plot(estTraj_Vanderpol(j, :), 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--'); hold on; grid on;
    legend('true', 'estimate');
end
title('vanderpol oscillator');


%% ========================LORENZ Model========================
disp("Lorenz example");
% Simulate True Trajectory & Measurement
for i_sim=2:sim.len
    trueTraj_Lorenz(:, i_sim) = Lorenz_process(trueTraj_Lorenz(:, i_sim-1), i_sim);
    [~, Z_Lorenz(:, i_sim)] = Identity_measurement(trueTraj_Lorenz(:, i_sim), 'R', sim.R_Lorenz);
end

% PF
for i_sim=2:sim.len
    % Propagation
    Ffcn_Lorenz = @(x)Lorenz_process(x, i_sim);
    pf_Lorenz_stadnard.Predict('TransitionFcn', Ffcn_Lorenz);
    
    % Measurement Update
    pf_Lorenz_stadnard.UpdateMeasurement(Z_Lorenz(:, i_sim));
    
    % Resample
    pf_Lorenz_stadnard.Resample('method', 'systematic', 'ratioNeff', sim.ratioNeff);

    estTraj_Lorenz(:, i_sim) = pf_Lorenz_stadnard.MMSE();
    disp([num2str(i_sim), ' th step done. ']);
end

figure;
for j=1:3
    subplot(3,1,j);
    plot(trueTraj_Lorenz(j, :), 'Color', colors(1,:), 'LineWidth', 1); hold on; grid on;
    plot(estTraj_Lorenz(j, :), 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--'); hold on; grid on;
    legend('true', 'estimate');
end
title('lorenz');

figure;
plot3(trueTraj_Lorenz(1, :), trueTraj_Lorenz(2, :), trueTraj_Lorenz(3, :), 'Color', colors(1,:), 'LineWidth', 1); hold on; grid on;
plot3(estTraj_Lorenz(1, :), estTraj_Lorenz(2, :), estTraj_Lorenz(3, :), 'Color', colors(2,:), 'LineWidth', 2, 'LineStyle', '--');
title('lorenz 3d');