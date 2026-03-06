% 清除环境
clear; close all; clc;

% 加载数据
load('PowerResult.mat');
Power = 40:-2:30; % 从40到30，步长为-2
ActiveAmp = 5:10;

% fontsize = 12;
figure('Position', [100, 100, 600, 500]);

% 创建网格
[A, P] = meshgrid(ActiveAmp, Power);

Prop = PropPowerNumSimuResult(end:-1:1, :);
Static = StaticPowerNumSimuResult(end:-1:1, :);
Passive = PassivePowerNumSimuResult(end:-1:1, :);

fontsize = 12;
% figure('Position', [100, 100, 600, 500]);

% 创建网格
% [P, A] = meshgrid(Power, ActiveAmp);

% === Passive ===
% === Static ===
baseColor = [200, 180, 240] / 255;  % 起始颜色
Zn = (Passive - min(Passive(:))) / (max(Passive(:)) - min(Passive(:)));
Zn = 1 - Zn;
R = baseColor(1) + (1 - baseColor(1)) * Zn;
G = baseColor(2) + (1 - baseColor(2)) * Zn;
B = baseColor(3) + (1 - baseColor(3)) * Zn;

% 构造 m×n×3 的 CData
CDataRGB = cat(3, R, G, B);

surf1 = surface(A, P, Passive, ...
    'EdgeColor', [200, 180, 240] / 255, ...               % 无曲面边线
    'FaceColor', [200, 180, 230] / 255, ...      % 自定义蓝色（非红色，但按你代码）
    'FaceAlpha', 'flat', ...
    'AlphaData', Passive, ...
    'LineWidth', 2.5);                       % 透明度由Z决定
hold on;
% === 生成反对角线：(i,j) -> (i+1, j-1) ===
x1 = A(1:end-1, 1:end-1);   % 起点 x
y1 = P(1:end-1, 1:end-1);   % 起点 y
z1 = Passive(1:end-1, 1:end-1);   % 起点 z

x2 = A(2:end, 2:end);   % 终点 x
y2 = P(2:end, 2:end);   % 终点 y
z2 = Passive(2:end, 2:end);   % 终点 z

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);

r1 = CDataRGB(1:end-1, 1:end-1, 1);
g1 = CDataRGB(1:end-1, 1:end-1, 2);
b1 = CDataRGB(1:end-1, 1:end-1, 3);

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);
r1 = r1(:); g1 = g1(:); b1 = b1(:);

% 用起点颜色画每条线（也可用中点色）
for k = 1:length(x1)
    plot3([x1(k), x2(k)], [y1(k), y2(k)], [z1(k), z2(k)], ...
        'Color', [200, 180, 240] / 255, ...
        'LineWidth', 2.5);
end

% === Static ===
baseColor = [169, 207, 131] / 255;  % 起始颜色
Zn = (Static - min(Static(:))) / (max(Static(:)) - min(Static(:)));
Zn = 1 - Zn;
R = baseColor(1) + (1 - baseColor(1)) * Zn;
G = baseColor(2) + (1 - baseColor(2)) * Zn;
B = baseColor(3) + (1 - baseColor(3)) * Zn;

% 构造 m×n×3 的 CData
CDataRGB = cat(3, R, G, B);

surf2 = surface(A, P, Static, ...
    'EdgeColor', 'interp', ...               % 无曲面边线
    'FaceColor', [180,210,150]/255, ...      % 自定义蓝色（非红色，但按你代码）
    'FaceAlpha', 'flat', ...
    'AlphaData', Static, ...
    'LineWidth', 2.5, ...
    'CData', CDataRGB, ...
    'CDataMapping', 'direct');                       % 透明度由Z决定
hold on;
% === 生成反对角线：(i,j) -> (i+1, j-1) ===
x1 = A(1:end-1, 1:end-1);   % 起点 x
y1 = P(1:end-1, 1:end-1);   % 起点 y
z1 = Static(1:end-1, 1:end-1);   % 起点 z

x2 = A(2:end, 2:end);   % 终点 x
y2 = P(2:end, 2:end);   % 终点 y
z2 = Static(2:end, 2:end);   % 终点 z

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);

r1 = CDataRGB(1:end-1, 1:end-1, 1);
g1 = CDataRGB(1:end-1, 1:end-1, 2);
b1 = CDataRGB(1:end-1, 1:end-1, 3);

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);
r1 = r1(:); g1 = g1(:); b1 = b1(:);

% 用起点颜色画每条线（也可用中点色）
for k = 1:length(x1)
    plot3([x1(k), x2(k)], [y1(k), y2(k)], [z1(k), z2(k)], ...
        'Color', [r1(k), g1(k), b1(k)], ...
        'LineWidth',2.5);
end

% === Prop ===
baseColor = [80, 191, 255] / 255;  % 起始颜色
Zn = (Prop - min(Prop(:))) / (max(Prop(:)) - min(Prop(:)));
Zn = 1 - Zn;
R = baseColor(1) + (1 - baseColor(1)) * Zn;
G = baseColor(2) + (1 - baseColor(2)) * Zn;
B = baseColor(3) + (1 - baseColor(3)) * Zn;

% 构造 m×n×3 的 CData
CDataRGB = cat(3, R, G, B);
% 设置透明度在 [0.2, 0.8] 范围内
alpha_min = 0.2;  % 最小透明度
alpha_max = 0.8;  % 最大透明度

alpha_data = alpha_min + (alpha_max - alpha_min) * ...
    (1 - (Prop - min(Prop(:))) / (max(Prop(:)) - min(Prop(:))));
surf3 = surface(A, P, Prop, ...
    'EdgeColor', 'interp', ...               % 无曲面边线
    'FaceColor', [80, 191, 255]/255, ...      % 自定义蓝色（非红色，但按你代码）
    'FaceAlpha', 'flat', ...
    'AlphaData', alpha_data, ...
    'LineWidth', 2.5, ...
    'CData', CDataRGB, ...
    'CDataMapping', 'direct');
hold on;
% === 生成反对角线：(i,j) -> (i+1, j-1) ===
x1 = A(1:end-1, 1:end-1);   % 起点 x
y1 = P(1:end-1, 1:end-1);   % 起点 y
z1 = Prop(1:end-1, 1:end-1);   % 起点 z

x2 = A(2:end, 2:end);   % 终点 x
y2 = P(2:end, 2:end);   % 终点 y
z2 = Prop(2:end, 2:end);   % 终点 z

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);

r1 = CDataRGB(1:end-1, 1:end-1, 1);
g1 = CDataRGB(1:end-1, 1:end-1, 2);
b1 = CDataRGB(1:end-1, 1:end-1, 3);

% 转为列向量
x1 = x1(:); y1 = y1(:); z1 = z1(:);
x2 = x2(:); y2 = y2(:); z2 = z2(:);
r1 = r1(:); g1 = g1(:); b1 = b1(:);

% 用起点颜色画每条线（也可用中点色）
for k = 1:length(x1)
    plot3([x1(k), x2(k)], [y1(k), y2(k)], [z1(k), z2(k)], ...
        'Color', [r1(k), g1(k), b1(k)], ...
        'LineWidth', 2.5);
end

% 更新 ColorData 为归一化值（double 类型）
% 设置标签
fontsize = 11;

% 设置标签并获取句柄
hX = xlabel('Maximum amplification factor, $\beta$', ...
    'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');

hY = ylabel('Maximum transmit power, $P_{\mathrm{t}}$ (dBm)', ...
    'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');

hZ = zlabel('Minimum SINR at UEs (dB)', ...
    'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');

% 设置视角
view(-45, 40);

% 设置字体和网格
set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman');
xlim([5, 10]);
ylim([30, 40]);
zlim([-15, 15]);
grid on;
% set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.8);
xticks(5:10);
yticks(30:2:40);
zticks(-15:5:15);
set(gcf, 'Color', 'white');
% 优化布局
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
disp(1)