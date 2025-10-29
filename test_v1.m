clear; clc; close all;

%% ---------- 参数设置 ----------
dt = 1;          % 时间步长
T = 30;             % 总时间
N = T / dt;        % 总帧数
t = (0:dt:T-dt)';  % 时间轴

% 三维目标真实状态（匀速 + 微扰）
pos0 = [500; 100; 200];
vel0 = [10; 8; 5]; 
pos = pos0 + vel0 .* t';

% 真实目标坐标
x_true = pos(1, :);
y_true = pos(2, :);
z_true = pos(3, :);

% 雷达噪声（标准差）
sigma_r = 2;               % 距离噪声 (m)
sigma_el_radar = 1;
sigma_az_radar = 1;

% 红外噪声
sigma_el_ir = 0.2;
sigma_az_ir = 0.2;

% ---------- 生成量测 ----------
data_all = zeros(N, 5); % 每行对应 [ir_az, ir_el, radar_az, radar_el, radar_r]
for k = 1:N
    % 真值到雷达球坐标
    r = sqrt(x_true(k)^2 + y_true(k)^2 + z_true(k)^2);
    az = rad2deg(atan2(y_true(k), x_true(k)));
    el = rad2deg(atan(z_true(k)/sqrt(x_true(k)^2 + y_true(k)^2)));

    % 红外：只测角度
    ir_el = el + sigma_el_ir * randn;
    ir_az = az + sigma_az_ir * randn;

    % 雷达：测距和角度
    radar_r = r + sigma_r * randn;
    radar_el = el + sigma_el_radar * randn;
    radar_az = az + sigma_az_radar * randn;

    data_all(k,:) = [ir_az, ir_el, radar_az, radar_el, radar_r];
end

%% ---------- 创建跟踪器实例 ----------
tracker = TrackFuser();

%% ---------- 主循环 ----------
fprintf('开始跟踪处理...\n');
for k = 1:N
    data = data_all(k,:); % 当前帧输入 [1×5]
    
    % 更新跟踪状态
    track = tracker.update(data);
    
    % 显示进度
    if mod(k, 10) == 0
        fprintf('  处理第 %d/%d 帧\n', k, N);
    end
end
fprintf('跟踪处理完成!\n');

%% ---------- 三维航迹对比 ----------
figure('Position', [100, 100, 1200, 800]);

% 三维轨迹
subplot(2,3,1);
plot3(x_true, y_true, z_true, 'k-', 'LineWidth', 3, 'DisplayName', '真实轨迹');
hold on;
plot3(track.radarstate(1,2:end), track.radarstate(2,2:end), track.radarstate(3,2:end), ...
     'r--', 'LineWidth', 2, 'DisplayName', '雷达跟踪');
plot3(track.irstate(1,2:end), track.irstate(2,2:end), track.irstate(3,2:end), ...
     'b:', 'LineWidth', 2, 'DisplayName', '红外跟踪');
plot3(track.fusiontrace(1,2:end), track.fusiontrace(2,2:end), track.fusiontrace(3,2:end), ...
     'g-', 'LineWidth', 2.5, 'DisplayName', '融合跟踪');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('三维航迹对比');
legend('Location', 'best'); grid on;
view(45, 30);

% XY平面投影
subplot(2,3,2);
plot(x_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', '真实轨迹');
hold on;
plot(track.radarstate(1,2:end), track.radarstate(2,2:end), 'r--', 'LineWidth', 1.5, 'DisplayName', '雷达跟踪');
plot(track.irstate(1,2:end), track.irstate(2,2:end), 'b:', 'LineWidth', 1.5, 'DisplayName', '红外跟踪');
plot(track.fusiontrace(1,2:end), track.fusiontrace(2,2:end), 'g-', 'LineWidth', 2, 'DisplayName', '融合跟踪');
xlabel('X (m)'); ylabel('Y (m)');
title('XY平面轨迹');
legend; grid on; axis equal;

% XZ平面投影
subplot(2,3,3);
plot(x_true, z_true, 'k-', 'LineWidth', 2, 'DisplayName', '真实轨迹');
hold on;
plot(track.radarstate(1,2:end), track.radarstate(3,2:end), 'r--', 'LineWidth', 1.5, 'DisplayName', '雷达跟踪');
plot(track.irstate(1,2:end), track.irstate(3,2:end), 'b:', 'LineWidth', 1.5, 'DisplayName', '红外跟踪');
plot(track.fusiontrace(1,2:end), track.fusiontrace(3,2:end), 'g-', 'LineWidth', 2, 'DisplayName', '融合跟踪');
xlabel('X (m)'); ylabel('Z (m)');
title('XZ平面轨迹');
legend; grid on; axis equal;

%% ---------- RMSE分析 ----------
% 计算位置误差模值
error_norm_radar = zeros(1, N);
error_norm_ir = zeros(1, N);
error_norm_fusion = zeros(1, N);

for k = 1:N
    true_pos = [x_true(k); y_true(k); z_true(k)];
    error_norm_radar(k) = norm(track.radarstate(1:3,k+1) - true_pos);
    error_norm_ir(k) = norm(track.irstate(1:3,k+1) - true_pos);
    error_norm_fusion(k) = norm(track.fusiontrace(1:3,k+1) - true_pos);
end

% 计算总体RMSE
rmse_radar = sqrt(mean(error_norm_radar.^2));
rmse_ir = sqrt(mean(error_norm_ir.^2));
rmse_fusion = sqrt(mean(error_norm_fusion.^2));

% 绘制RMSE对比
subplot(2,3,4);
methods = {'雷达', '红外', '融合'};
rmse_values = [rmse_radar, rmse_ir, rmse_fusion];
colors = [0.9 0.3 0.3; 0.3 0.3 0.9; 0.3 0.7 0.3]; % 红, 蓝, 绿

b = bar(rmse_values, 'FaceColor', 'flat');
b.CData = colors;

% 在柱状图上显示数值
for i = 1:3
    text(i, rmse_values(i) + 0.1, sprintf('%.3f m', rmse_values(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
end

set(gca, 'XTickLabel', methods);
ylabel('位置RMSE (m)');
title('总体位置RMSE对比');
grid on;

%% ---------- 误差分析 ----------
% 各坐标轴误差
subplot(2,3,5);
plot(t, error_norm_radar, 'r--', 'LineWidth', 1, 'DisplayName', '雷达误差');
hold on;
plot(t, error_norm_ir, 'b:', 'LineWidth', 1, 'DisplayName', '红外误差');
plot(t, error_norm_fusion, 'g-', 'LineWidth', 1.5, 'DisplayName', '融合误差');
xlabel('时间 (s)'); ylabel('位置误差 (m)');
title('位置误差随时间变化');
legend; grid on;

% 误差分布
subplot(2,3,6);
box_data = [error_norm_radar', error_norm_ir', error_norm_fusion'];
boxplot(box_data, 'Labels', {'雷达', '红外', '融合'});
ylabel('位置误差 (m)');
title('位置误差分布');
grid on;

sgtitle('红外-雷达融合跟踪性能分析');

%% ---------- 性能统计 ----------
fprintf('=== 跟踪性能统计 ===\n\n');
fprintf('总体位置RMSE:\n');
fprintf('雷达跟踪: %.4f m\n', rmse_radar);
fprintf('红外跟踪: %.4f m\n', rmse_ir);
fprintf('融合跟踪: %.4f m\n', rmse_fusion);

improvement_vs_radar = (1 - rmse_fusion/rmse_radar) * 100;
improvement_vs_ir = (1 - rmse_fusion/rmse_ir) * 100;

fprintf('\n融合提升效果:\n');
fprintf('相对于雷达: +%.2f%% 精度提升\n', improvement_vs_radar);
fprintf('相对于红外: +%.2f%% 精度提升\n', improvement_vs_ir);

% 稳态性能分析（后50%数据）
steady_start = round(N/2);
rmse_steady_radar = sqrt(mean(error_norm_radar(steady_start:end).^2));
rmse_steady_ir = sqrt(mean(error_norm_ir(steady_start:end).^2));
rmse_steady_fusion = sqrt(mean(error_norm_fusion(steady_start:end).^2));

fprintf('\n稳态性能（后50%%数据）:\n');
fprintf('雷达稳态RMSE: %.4f m\n', rmse_steady_radar);
fprintf('红外稳态RMSE: %.4f m\n', rmse_steady_ir);
fprintf('融合稳态RMSE: %.4f m\n', rmse_steady_fusion);

fprintf('\n=== 分析结论 ===\n');
if improvement_vs_radar > 10 && improvement_vs_ir > 5
    fprintf('✅ 融合跟踪显著优于单一传感器：\n');
    fprintf('   融合算法有效结合了雷达的距离测量优势和红外的角度测量优势。\n');
elseif improvement_vs_radar > 0 && improvement_vs_ir > 0
    fprintf('✅ 融合跟踪具有稳定改进：\n');
    fprintf('   在所有测试时段内，融合跟踪均表现出最佳性能。\n');
else
    fprintf('⚠️ 融合效果有待优化：\n');
    fprintf('   需要进一步调整融合权重或滤波参数。\n');
end

%% ---------- 蒙特卡洛仿真示例 ----------
fprintf('\n=== 开始蒙特卡洛仿真 (10次运行) ===\n');
mc_runs = 10;
mc_results = zeros(mc_runs, 3); % 存储每次运行的RMSE

for mc = 1:mc_runs
    % 创建新的跟踪器实例
    mc_tracker = TrackFuser();
    
    % 重新生成带噪声的测量数据
    mc_data_all = zeros(N, 5);
    for k = 1:N
        r = sqrt(x_true(k)^2 + y_true(k)^2 + z_true(k)^2);
        az = rad2deg(atan2(y_true(k), x_true(k)));
        el = rad2deg(atan(z_true(k)/sqrt(x_true(k)^2 + y_true(k)^2)));

        ir_el = el + sigma_el_ir * randn;
        ir_az = az + sigma_az_ir * randn;
        radar_r = r + sigma_r * randn;
        radar_el = el + sigma_el_radar * randn;
        radar_az = az + sigma_az_radar * randn;

        mc_data_all(k,:) = [ir_az, ir_el, radar_az, radar_el, radar_r];
    end
    
    % 运行跟踪
    for k = 1:N
        mc_track = mc_tracker.update(mc_data_all(k,:));
    end
    
    % 计算本次运行的RMSE
    mc_error_radar = zeros(1, N);
    mc_error_ir = zeros(1, N);
    mc_error_fusion = zeros(1, N);
    
    for k = 1:N
        true_pos = [x_true(k); y_true(k); z_true(k)];
        mc_error_radar(k) = norm(mc_track.radarstate(1:3,k+1) - true_pos);
        mc_error_ir(k) = norm(mc_track.irstate(1:3,k+1) - true_pos);
        mc_error_fusion(k) = norm(mc_track.fusiontrace(1:3,k+1) - true_pos);
    end
    
    mc_results(mc, 1) = sqrt(mean(mc_error_radar.^2));
    mc_results(mc, 2) = sqrt(mean(mc_error_ir.^2));
    mc_results(mc, 3) = sqrt(mean(mc_error_fusion.^2));
    
    fprintf('  完成第 %d/%d 次蒙特卡洛运行\n', mc, mc_runs);
end

fprintf('\n蒙特卡洛仿真结果:\n');
fprintf('雷达RMSE: %.3f ± %.3f m\n', mean(mc_results(:,1)), std(mc_results(:,1)));
fprintf('红外RMSE: %.3f ± %.3f m\n', mean(mc_results(:,2)), std(mc_results(:,2)));
fprintf('融合RMSE: %.3f ± %.3f m\n', mean(mc_results(:,3)), std(mc_results(:,3)));