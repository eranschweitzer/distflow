%%%% view statistics of the alpha estimation
clear variables
close all
%% load data
load('alphaest_opt.mat')

%% plot mean and median
figure;
plot(feeder_sizes, stats.mean, 'k', 'linewidth', 2)
hold on;
plot(feeder_sizes, stats.mean + stats.std, '--k', 'linewidth', 1.5)
plot(feeder_sizes, stats.mean - stats.std, '--k', 'linewidth', 1.5)
xlabel('# of nodes')
ylabel('\alpha')
title('mean \pm std')

figure;
plot(feeder_sizes, stats.median, 'k', 'linewidth', 2)
title('median')
%% plot histograms
idx = [1, round(length(feeder_sizes)/2), length(feeder_sizes)];
figure;
for k = 1:3
    subplot(3,1,k)
    histogram(alpha{idx(k)},'Normalization', 'probability')
    xlabel('\alpha')
    ylabel('Probability')
    titlestr = sprintf('Feeders with %d Nodes', feeder_sizes(idx(k)));
    title(titlestr)
end