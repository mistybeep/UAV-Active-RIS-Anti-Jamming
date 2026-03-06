% 清除环境
clear; close all; clc;

% 加载数据
load('aaa.mat');
% load("Prop_StaticRISNumRes.mat")

% 提取数据
Prop = PropRISNumSimuRes(end,:);
MM = MMRISNumSimuRes(end,:);
SDR = SDRRISNumSimuRes(end,:);
Static = PropStaticRISNumSimuRes(end,:);
Passive = PassiveRISNumSimuRes(end,:);
load('ErrorNumberRes.mat');
Prop(9:11) = PropRISNumSimuRes(end,9:11);
MM(9:11) = MMRISNumSimuRes(end,9:11);
SDR(9:11) = SDRRISNumSimuRes(end,9:11);
Static(9:11) = PropStaticRISNumSimuRes(end,9:11);
Passive(9:11) = PassiveRISNumSimuRes(end,9:11);
% 创建x轴数据
iter = 0:2:20;

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
     'MarkerSize', 8, 'LineWidth', 1.2, ...
     'DisplayName', 'Proposed Method');

% 然后绘制其他线
h(2) = plot(iter, MM, '-', 'Color', C(2,:), ...
     'Marker', 'd', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(2,:), ...
     'MarkerSize', 7, 'LineWidth', 1.2, ...
     'DisplayName', 'MM (Benchmark)');

h(3) = plot(iter, SDR, '-', 'Color', C(3,:), ...
     'Marker', 'o', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(3,:), ...
     'MarkerSize', 7, 'LineWidth', 1.2, ...
     'DisplayName', 'SDR (Benchmark)');

h(4) = plot(iter, Static, '-', 'Color', C(4,:), ...
     'Marker', '^', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(4,:), ...
     'MarkerSize', 7, 'LineWidth', 1.2, ...
     'DisplayName', 'Static Active RIS');

h(5) = plot(iter, Passive, '-', 'Color', C(5,:), ...
     'Marker', 'v', ...
     'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(5,:), ...
     'MarkerSize', 7, 'LineWidth', 1.2, ...
     'DisplayName', 'Passive RIS');

% 设置坐标轴标签和其他属性
xlabel('Jamming Position Uncertainty, $\epsilon$ (m)', 'Interpreter', 'latex', ...
       'FontSize', fontsize, 'FontName', 'Times New Roman');
ylabel('Minimum SINR at UEs (dB)', 'FontSize', fontsize, ...
       'FontName', 'Times New Roman');
xlim([0, 20]);
ylim([-20, 0]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
xticks(0:2:20);
yticks(-20:5:0);
set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman');

% 手动调整显示顺序：将Prop线移到最前面
uistack(h(1), 'top');

% 创建图例（按原始顺序）
legend(h, 'Proposed Method', 'MM (Benchmark)', 'SDR (Benchmark)', 'Static Active RIS', 'Passive RIS', ...
       'Location', 'southeast', 'Box', 'on', 'EdgeColor', 'black', ...
       'FontSize', fontsize-1.5, 'NumColumns', 1);
set(gca, 'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.5); 
set(gcf, 'Color', 'white');
set(gcf, 'Renderer', 'painters'); 
% 2. 导出图形
exportgraphics(gcf, 'Error.pdf', ...
    'ContentType', 'vector');