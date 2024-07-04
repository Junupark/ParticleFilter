sim_len = 1000;
trueX = 100*ones(3,1);
trueTraj = zeros(length(trueX), sim_len);
trueTraj(:,1) = trueX;

for i = 2:sim_len
    trueTraj(:,i) = Lorenz_process(trueTraj(:,i-1), i);
end
figure; plot3(trueTraj(1,:), trueTraj(2,:), trueTraj(3,:)); grid on;
% figure; plot(trueTraj(1,:), trueTraj(2,:), 'ko'); grid on;
% figure; plot(1:sim_len, trueTraj(1,:)); grid on;


[t,x] = ode45(@(t,x)[10*(x(2) - x(1)); x(1)*(28 - x(3)) - x(2); x(1)*x(2) - 8/3*x(3)], [0, 0.05*sim_len], 100*ones(3,1));
hold on; plot3(x(:,1), x(:,2), x(:,3), 'r-');
