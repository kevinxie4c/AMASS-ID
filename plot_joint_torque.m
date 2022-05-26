forces = load('forces.txt');
t = 1:1/1200:100;
len = 400;
plot(t(1:len),forces(1:len,10),t(1:len),forces(1:len,22), 'LineWidth', 3);
xlabel('time (s)')
ylabel('torque (N \cdot m)')
legend('Left knee', 'Right knee')
title('Torques on knees')
