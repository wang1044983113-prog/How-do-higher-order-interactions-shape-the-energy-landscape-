% 简化版：复现论文图2的核心特征
% 展示从非扭曲态通过噪声转变到扭曲态的过程
clc, clear, close all;

% 参数设置
n = 83;
r = 2;
sigma = 1;
sigma_delta = 3.5;
D = 0.4;

T = 400;
dt = 0.01;
steps = T/dt;
noise_stop_time = 300;
noise_stop_step = noise_stop_time/dt;

%% 搜索非扭曲态初始条件
fprintf('正在搜索非扭曲态初始条件...\n');
found = false;

for attempt = 1:10
    test_init = 2*pi*rand(n, 1) - pi;
    test_theta = test_init;
    
    % 演化50秒看是否稳定到非扭曲态
    for t = 1:1000
        dtheta = compute_interactions(test_theta, n, r, sigma, sigma_delta);
        test_theta = test_theta + dtheta * dt;
        test_theta = mod(test_theta + pi, 2*pi) - pi;
    end
    
    % 检查是否为非扭曲态
    local_order = calculate_local_order(test_theta, n, r);
    p_twisted = sum(local_order >= 0.85) / n;
    
    if p_twisted < 0.95
        initial_phases = test_init;
        found = true;
        fprintf('找到非扭曲态！P_twisted = %.4f (attempt %d)\n', p_twisted, attempt);
        break;
    end
end

if ~found
    fprintf('未找到合适初始条件，使用随机初始条件\n');
    initial_phases = 2*pi*rand(n, 1) - pi;
end

%% 两个系统的并行模拟
fprintf('开始模拟...\n');

% 初始化
theta_noise = initial_phases;
theta_clean = initial_phases;

% 历史记录
p_noise_hist = zeros(1, steps);
p_clean_hist = zeros(1, steps);
theta_noise_hist = zeros(n, steps);

% 初始值
local_order = calculate_local_order(theta_noise, n, r);
p_noise_hist(1) = sum(local_order >= 0.85) / n;
p_clean_hist(1) = p_noise_hist(1);
theta_noise_hist(:, 1) = theta_noise;

% 主循环
for t = 2:steps
    % 有噪声系统
    dtheta = compute_interactions(theta_noise, n, r, sigma, sigma_delta);
    theta_noise = theta_noise + dtheta * dt;
    
    % 只在前300秒加噪声
    if (t-1)*dt < noise_stop_time
        theta_noise = theta_noise + sqrt(2*D*dt) * randn(n, 1);
    end
    
    theta_noise = mod(theta_noise + pi, 2*pi) - pi;
    theta_noise_hist(:, t) = theta_noise;
    
    local_order = calculate_local_order(theta_noise, n, r);
    p_noise_hist(t) = sum(local_order >= 0.85) / n;
    
    % 无噪声系统
    dtheta = compute_interactions(theta_clean, n, r, sigma, sigma_delta);
    theta_clean = theta_clean + dtheta * dt;
    theta_clean = mod(theta_clean + pi, 2*pi) - pi;
    
    local_order = calculate_local_order(theta_clean, n, r);
    p_clean_hist(t) = sum(local_order >= 0.85) / n;
    
    if mod(t, 5000) == 0
        fprintf('进度: %.1f%%\n', 100*t/steps);
    end
end

%% 噪声移除后的5秒检查
fprintf('执行5秒噪声移除检查...\n');
theta_check = theta_noise_hist(:, noise_stop_step);

for t = 1:10000  
    dtheta = compute_interactions(theta_check, n, r, sigma, sigma_delta);
    theta_check = theta_check + dtheta * dt;
    theta_check = mod(theta_check + pi, 2*pi) - pi;
end

local_order_check = calculate_local_order(theta_check, n, r);
p_check = sum(local_order_check >= 0.85) / n;
fprintf('移除噪声5秒后: P_twisted = %.4f\n', p_check);

%% 绘图
figure('Position', [100, 100, 1400, 500], 'Color', 'white');

% 子图1: P_twisted时间演化
subplot(1, 3, 1);
time = (0:steps-1)*dt;
plot(time, p_noise_hist, 'r-', 'LineWidth', 2);
hold on;
plot(time, p_clean_hist, 'b-', 'LineWidth', 2);
xline(noise_stop_time, 'k--', 'LineWidth', 1.5);
text(noise_stop_time+10, 0.5, 'Noise removed', 'FontSize', 11);
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('P_{twisted}', 'FontSize', 12, 'Interpreter', 'tex');
title('(a) Temporal Evolution', 'FontSize', 13, 'FontWeight', 'bold');
legend('With Noise', 'Without Noise', 'Location', 'best');
grid on;
ylim([0 1.05]);
xlim([0 T]);
box on;

% 子图2: t=300s的相位分布
subplot(1, 3, 2);
scatter(1:n, theta_noise_hist(:, noise_stop_step), 60, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Oscillator Index', 'FontSize', 12);
ylabel('Phase', 'FontSize', 12);
title('(b) Phase at t=300s', 'FontSize', 13, 'FontWeight', 'bold');
ylim([-pi pi]);
yticks(-pi:pi/2:pi);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
grid on;
box on;

% 子图3: 移除噪声5秒后的相位分布
subplot(1, 3, 3);
scatter(1:n, theta_check, 60, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Oscillator Index', 'FontSize', 12);
ylabel('Phase', 'FontSize', 12);
title('(c) After 5s noise-free', 'FontSize', 13, 'FontWeight', 'bold');
ylim([-pi pi]);
yticks(-pi:pi/2:pi);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
grid on;
box on;

sgtitle(sprintf('Figure 2 Reproduction: σ=%.1f, σ_Δ=%.1f, D=%.1f', ...
    sigma, sigma_delta, D), 'FontSize', 14, 'FontWeight', 'bold');

%% 辅助函数
function dtheta = compute_interactions(theta, n, r, sigma, sigma_delta)
    dtheta = zeros(n, 1);
    for i = 1:n
        % 成对相互作用
        pair_sum = 0;
        for j = i-r:i+r
            j_idx = mod(j-1, n) + 1;
            if j_idx ~= i
                pair_sum = pair_sum + sin(theta(j_idx) - theta(i));
            end
        end
        
        % 三元相互作用
        triad_sum = 0;
        for j = i-r:i+r
            for k = i-r:i+r
                j_idx = mod(j-1, n) + 1;
                k_idx = mod(k-1, n) + 1;
                if j_idx ~= i && k_idx ~= i && j_idx ~= k_idx
                    triad_sum = triad_sum + sin(theta(j_idx) + theta(k_idx) - 2*theta(i));
                end
            end
        end
        
        dtheta(i) = sigma/(2*r) * pair_sum + sigma_delta/(2*r*(2*r-1)) * triad_sum;
    end
end

function local_order = calculate_local_order(theta, n, r)
    local_order = zeros(n, 1);
    for j = 1:n
        local_sum = 0;
        for k = j-r:j+r
            k_idx = mod(k-1, n) + 1;
            local_sum = local_sum + exp(1i * theta(k_idx));
        end
        local_order(j) = abs(local_sum) / (2*r+1);
    end
end