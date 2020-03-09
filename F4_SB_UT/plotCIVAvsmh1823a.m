% Plot PODs

clear; clc;

% Load CIVA POD.
podCIVA = load('podCIVA');

% Get mh1823a POD.
pod1823 = f4sb_pod_glm(...
    'logX', 0, ...
    'logY', 0, ...
    'ahatDecFactor', 0.2, ...
    'dispPlot', 0);

% Plot
h(1) = plot(podCIVA.aCIVA*25.4, podCIVA.pod50CIVA, '-r');hold on
h(2) = plot(podCIVA.aCIVA*25.4, podCIVA.pod95CIVA, '--r');
h(3) = plot(pod1823.aPOD50, pod1823.pod, '-b');
h(4) = plot(pod1823.aPOD95, pod1823.pod, '--b'); hold off

% Format
grid on
axis([0.2, 1.4, 0, 1]);
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');

% legend
leg = legend(h([1, 3]), {'CIVA', 'mh1823a'});
set(leg, 'Location', 'northwest');
set(leg, 'Color', 'none');
set(leg, 'FontSize', 14);
set(leg, 'FontWeight', 'bold');

% Label
xl = xlabel('Crack Length (mm)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Probability of Detection');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');