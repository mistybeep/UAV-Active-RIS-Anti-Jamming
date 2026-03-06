% 清除环境
clear; close all; clc;

% 加载数据
load('RISNumberRes.mat');
% load("Prop_StaticRISNumRes.mat")

% 提取数据
Prop = PropRISNumSimuRes(end,:);
MM = MMRISNumSimuRes(end,:);
SDR = SDRRISNumSimuRes(end,:);
Static = PropStaticRISNumSimuRes(end,:);
Passive = PassiveRISNumSimuRes(end,:);

% 创建x轴数据
iter = [2:9].^2;

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
xlabel('Number of Active RIS elements, $N$', 'Interpreter', 'latex', ...
       'FontSize', fontsize, 'FontName', 'Times New Roman');
ylabel('Minimum SINR at UEs (dB)', 'FontSize', fontsize, ...
       'FontName', 'Times New Roman');
xlim([0, 90]);
ylim([-30, 10]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
xticks(0:10:100);
yticks(-30:5:10);
set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman');

% 手动调整显示顺序：将Prop线移到最前面
uistack(h(1), 'top');

% 创建图例（按原始顺序）
legend(h, 'Proposed Method', 'MM (Benchmark)', 'SDR (Benchmark)', 'Static Active RIS', 'Passive RIS', ...
       'Location', 'southeast', 'Box', 'on', 'EdgeColor', 'black', ...
       'FontSize', fontsize-1.5, 'NumColumns', 1);
set(gca, 'TickDir', 'in', 'Box', 'on', 'LineWidth', 0.5); 
set(gcf, 'Color', 'white');

axesPosition = [0.33 0.17 0.24 0.21]; 
axes('Position', axesPosition);
hold on;

% 2. 在小图中重新绘制数据
% 这里稍微减小了线宽和标记点大小，让小图看起来更精致

plot(iter, MM, '-', 'Color', C(2,:), 'Marker', 'd', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(2,:), 'MarkerSize', 6, 'LineWidth', 1.1);
plot(iter, SDR, '-', 'Color', C(3,:), 'Marker', 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(3,:), 'MarkerSize', 6, 'LineWidth', 1.1);
plot(iter, Static, '-', 'Color', C(4,:), 'Marker', '^', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(4,:), 'MarkerSize', 6, 'LineWidth', 1.1);
plot(iter, Passive, '-', 'Color', C(5,:), 'Marker', 'v', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(5,:), 'MarkerSize', 6, 'LineWidth', 1.1);
plot(iter, Prop, '-s', 'Color', C(1,:), 'MarkerFaceColor', 'white', 'MarkerEdgeColor', C(1,:), 'MarkerSize', 7, 'LineWidth', 1.1);
% 3. 设置你想要放大的数据区间 (X轴和Y轴范围)
% 【注意】你需要根据你的实际数据修改这里的区间！
zoom_x = [63, 82]; % 假设你想放大 N 在 60 到 85 之间的部分
zoom_y = [5, 8]; % 如果需要，也可以限制 Y 轴。如果不写，MATLAB会自动适应 X 轴区间的 Y 值
xlim(zoom_x);
ylim(zoom_y); 

% 4. 设置小图的样式
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
set(gca, 'FontSize', 9, 'FontName', 'Times New Roman'); % 字体比主图小一点
box on; % 给小图加个边框