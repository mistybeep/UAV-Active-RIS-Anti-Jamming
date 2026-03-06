% 清除环境
clear; close all; clc;

% 加载数据
load('JamPowerRes.mat');
% load('Prop_StaticJamPowNumRes.mat');
% 提取数据并转换为行向量
Prop = PropJamPowSimuRes(end, :)';
MM = MMJamPowSimuRes(end, :)';
SDR = SDRJamPowSimuRes(end, :)';
Static = PropStaticJamPowSimuRes(end, :)';
Passive = PassiveJamPowSimuRes(end, :)';
Pow = [1, 3];
% 提取前2个元素
PropData = Prop(Pow);
MMData = MM(Pow);
SDRData = SDR(Pow);
StaticData = Static(Pow);
PassiveData = Passive(Pow);
% 设置图形参数
fontsize = 12;
figure('Position', [100, 100, 500, 400]);

% 设置x轴标签
JamPower = {'$\mathrm{SJNR = -20~dB}$', '$\mathrm{SJNR = 0~dB}$'};
x = 0:length(JamPower)-1;  % [0, 1]
width = 0.1;
X1 = x;
X2 = x + width + 0.02;
X3 = X2 + width + 0.02;
X4 = X3 + width + 0.02;
X5 = X4 + width + 0.02;

center_positions = X3;

% 定义颜色
colors = [120, 148, 84;    % #789454 - 橄榄绿
          255, 214, 54;    % #FFD636 - 金黄色
          94, 135, 186;    % #5E87BA - 灰蓝色
          212, 143, 128;
          154, 123, 176] / 255; % #D48F80 - 珊瑚粉

% 推荐的第五个颜色
% #9A7BB0 - 淡紫色 (154, 123, 176)
fifth_color = [154, 123, 176] / 255;

PassiveData = PassiveData + 4;

% 绘制柱状图
bar(X1, PropData, width, 'FaceColor', colors(1,:), 'DisplayName', 'Proposed Method');
hold on;
bar(X2, MMData, width, 'FaceColor', colors(2,:), 'DisplayName', 'MM (Benchmark)');
bar(X3, SDRData, width, 'FaceColor', colors(3,:), 'DisplayName', 'SDR (Benchmark)');
bar(X4, StaticData, width, 'FaceColor', colors(4,:), 'DisplayName', 'Static Active RIS');
bar(X5, PassiveData, width, 'FaceColor', colors(5,:), 'DisplayName', 'Passive RIS');

% 设置坐标轴标签
xlabel('SJNR (dB)', 'Interpreter', 'latex', ...
       'FontSize', fontsize, 'FontName', 'Times New Roman');
ylabel('Minimum SINR at UEs (dB)', 'FontSize', fontsize, ...
       'FontName', 'Times New Roman');

% 设置x轴刻度
xticks(center_positions);
xticklabels(JamPower);
% yticks(-10:2:10);
set(gca, 'TickLabelInterpreter', 'latex');

% 设置字体
set(gca, 'FontSize', fontsize, 'FontName', 'Times New Roman');
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6, 'GridColor', [0.5 0.5 0.5]);
% 设置图例
legend('Location', 'northeast', ...
       'Box', 'on', ...
       'EdgeColor', 'black', ...
       'FontSize', fontsize-2,  ...
       'NumColumns', 1, ...
       'Box', 'on');

% 优化布局
set(gcf, 'Color', 'white');
set(gca, 'TickDir', 'out');
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
t1=text(0, 0.23,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t1,'rotation',90);
t2=text(1, 0.23,'//','sc','BackgroundColor','w','margin',eps, 'fontsize',13);
set(t2,'rotation',90);
% 调整坐标轴范围（如果需要）
xlim([min(X1)-2*width, max(X5)+2*width]);
ylim([-5, 7]);
ytick=-5:2:7;
set(gca,'ytick',ytick);
ytick(ytick<-2-eps)=ytick(ytick<-2-eps)-4;
for i=1:length(ytick)
   yticklabel{i}=sprintf('%d',ytick(i));
end
set(gca,'yTickLabel', yticklabel, 'FontSize', 12, 'FontName', 'Times New Roman'); %修改坐标名称、字体
% 显示图形
hold off;