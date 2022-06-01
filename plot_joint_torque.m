forces = load('output/forces.txt');
error = load('output/errors.txt');
t = 1:1/1200:100;
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
