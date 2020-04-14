function [] = olhsCIVADesign()
% OLHSCIVADESIGN is a function to get experimental conditions for CIVA
% simulation using optimal latin hypercube sampling (OLHS).
%
% Inputs: none
% Outputs: none
%
% Revision history
% 022516 LDY Code is modified from 'olhsSIFDesign.m'.

% Constants
t = 4.08;    % Tang thickness, [mm]

% -------- Optimal Latin Hypercube Sampling (50 Samples) --------
x = lhsdesign(50, 6, 'criterion', 'maximin');
freq = [0.8, 2];    % Frequency, [MHz]
r1 = [0.1, t/2];    % Half width, [mm]
r2 = [0.1, t];    % Depth, [mm]
sPath = [0, t];    % Scan path, [mm]
liftOff = [0, 1];    % Lift-off, [mm]
cond = [0.4, 1];    % Conductivity, [%IACS]
cond = cond*58.001*0.01;    % Conductivity (100%IACS = 58.001MS/m), [MS/m]

x(:, 1) = x(:, 1)*(freq(2) - freq(1)) + freq(1);
x(:, 2) = x(:, 2)*(r1(2) - r1(1)) + r1(1);
x(:, 3) = x(:, 3)*(r2(2) - r2(1)) + r2(1);
x(:, 4) = x(:, 4)*(sPath(2) - sPath(1)) + sPath(1);
x(:, 5) = x(:, 5)*(liftOff(2) - liftOff(1)) + liftOff(1);
x(:, 6) = x(:, 6)*(cond(2) - cond(1)) + cond(1);

% Write csv file
headers = {'freq_MHz', 'r1_mm', 'r2_mm', 'sPath_mm', 'liftOff_mm', 'cond_MS'};
data = x;
csvwrite_with_headers([date, '_olhs_civa.csv'], data, headers);

figure(1)
plot(x(:, 2), x(:, 3), 'ro', 'MarkerFaceColor', 'b')
xlabel('r1 (mm)')
ylabel('r2 (mm)')
title('OLHS (50 Samples)')
set(gcf, 'color', 'w')


end