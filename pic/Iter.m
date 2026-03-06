% 清除环境
clear; close all; clc;

% 加载数据
load('RISNumberRes.mat');
% load("Prop_StaticRISNumRes.mat")

% 提取数据
Prop = PropRISNumSimuRes(:,3);
MM = MMRISNumSimuRes(:,3);
SDR = SDRRISNumSimuRes(:,3);
Static = PropStaticRISNumSimuRes(:,3);
Passive = PassiveRISNumSimuRes(:,3);
% Prop2 = PropRISNumSimuRes(:,6);
% MM2 = MMRISNumSimuRes(:,6);
% SDR2 = SDRRISNumSimuRes(:,6);
load('IterNumberRes.mat');
Prop2 = PropRISNumSimuRes(:,3);
MM2 = MMRISNumSimuRes(:,3);
SDR2 = SDRRISNumSimuRes(:,3);
% 创建x轴数据
iter = 0:10;

% 设置图形参数
fontsize = 12;
figure('Position', [100, 100, 500, 400]);

hold on;
C = [
    1.0000    0.0000    0.0000    % 黑色 (No prediction)
    1.0000    0.5000    0.0000    % 橙色 (AR)
    0.0000    0.0000    1.0000    % 蓝色 (PVEC)
    1.0000    0.0000    1.0000    % 品红 (STEM kernel learning)
    0.0000    0.0000    0.0000    % 红色 (GEM kernel learning)
];
% 保存所有图形对象的句柄
h = [];

% 先绘制Prop线（显示在最上面）
h(1) = plot(iter, Prop, '-s', 'Color', C(1,:), ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(1,:), ...
     'MarkerSize', 8, 'LineWidth', 1.1, ...
     'DisplayName', 'Proposed Method');

% 然后绘制其他线
h(2) = plot(iter, MM, '-', 'Color', C(2,:), ...
     'Marker', 'd', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(2,:), ...
     'MarkerSize', 7, 'LineWidth', 1.1, ...
     'DisplayName', 'MM (Benchmark)');

h(3) = plot(iter, SDR, '-', 'Color', C(3,:), ...
     'Marker', 'o', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(3,:), ...
     'MarkerSize', 7, 'LineWidth', 1.1, ...
     'DisplayName', 'SDR (Benchmark)');

% h(4) = plot(iter, Prop2, '-s', 'Color', [1, 0, 0], ...
%      'MarkerFaceColor', 'white', 'MarkerEdgeColor', [1, 0, 0], ...
%      'MarkerSize', 8, 'LineWidth', 1.2, ...
%      'DisplayName', 'Proposed Method');
% 
% % 然后绘制其他线
% h(5) = plot(iter, MM2, '-', 'Color', [1, 0, 1], ...
%      'Marker', 'd', ...
%      'MarkerFaceColor', 'white', 'MarkerEdgeColor', [1, 0, 1], ...
%      'MarkerSize', 7, 'LineWidth', 1.2, ...
%      'DisplayName', 'MM (Benchmark)');
% 
% h(6) = plot(iter, SDR2, '-', 'Color', [0, 0, 1], ...
%      'Marker', 'o', ...
%      'MarkerFaceColor', 'white', 'MarkerEdgeColor', [0, 0, 1], ...
%      'MarkerSize', 7, 'LineWidth', 1.2, ...
%      'DisplayName', 'SDR (Benchmark)');


h(4) = plot(iter, Static, '-', 'Color', C(4,:), ...
     'Marker', '^', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(4,:), ...
     'MarkerSize', 7, 'LineWidth', 1.1, ...
     'DisplayName', 'Static Active RIS');
% 
% h(5) = plot(iter, Passive, '-', 'Color', [0, 0, 0], ...
%      'Marker', 'v', ...
%      'MarkerFaceColor', 'white', 'MarkerEdgeColor', [0, 0, 0], ...
%      'MarkerSize', 7, 'LineWidth', 1.2, ...
%      'DisplayName', 'Passive RIS');

% 设置坐标轴标签和其他属性
xlabel('Number of Iterations, $I$', 'Interpreter', 'latex', ...
       'FontSize', fontsize, 'FontName', 'Times New Roman');
ylabel('Minimum SINR at UEs (dB)', 'FontSize', fontsize, ...
       'FontName', 'Times New Roman');
xlim([0, 10]);
ylim([-25, 0]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
xticks(0:1:10);
yticks(-25:5:0);
set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman');

% 手动调整显示顺序：将Prop线移到最前面
uistack(h(1), 'top');

% 创建图例（按原始顺序）
legend(h, 'Proposed Method', 'MM (Benchmark)', 'SDR (Benchmark)', 'Static Active RIS', ...
       'Location', 'southeast', 'Box', 'on', 'EdgeColor', 'black', ...
       'FontSize', fontsize-1.5, 'NumColumns', 1);

set(gcf, 'Color', 'white');
% 优化布局
% set(gca, 'TickDir', 'in');
set(gca, 'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.5); 
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
% 保存图形
set(gcf, 'Renderer', 'painters'); 
% 2. 导出图形
exportgraphics(gcf, 'Iter.pdf', ...
    'ContentType', 'vector'); % 删掉 Padding 这一行，默认就是 tight