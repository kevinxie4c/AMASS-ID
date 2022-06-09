forces = load('output/forces.txt');
error = load('output/errors.txt');
% forces = load('jump-output/forces.txt');
% error = load('jump-output/errors.txt');
t = 0:1/120:100;
len = 400;

yyaxis left
plot(t(1:len),forces(1:len,10), 'r-', 'LineWidth', 3);
hold on
plot(t(1:len),forces(1:len,22), 'g-', 'LineWidth', 3);
xlabel('time (s)')
ylabel('torque (N \cdot m)')

yyaxis right
plot(t(1:len), error(1:len), 'b-', 'LineWidth', 3);
ylabel('error');
ylim([0, 600]);

legend('Left knee', 'Right knee', 'error')
title('Torques on knees')

forces2 = load('walking_base_tau_filtered.txt');
tt = 0:1/600:(760-1)/600;
n1 = 152;
n2 = 760;
figure;
plot(t(1:n1),forces(105:104+n1,10), 'LineWidth', 3);
hold on;
plot(tt(1:n2),forces2(1:n2,20), 'LineWidth', 3);

