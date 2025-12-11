clc, clear, close all;

% 记录开始时间
tic;

% 参数设置
n = 83;                  % 振荡器数量
r = 2;                   % 耦合范围
sigma = 1;               % 成对耦合强度（已通过时间缩放设为1）
sigma_delta = 3;         % 三元耦合强度，固定值
D = 0.4;                 % 较小的噪声强度，保持稳定态

% 模拟参数
T = 30000;                 % 总模拟时间
dt = 0.1;                % 时间步长
steps = T/dt;

% 定义要存储的时间点
time_points = 1:1:steps;   % 每10步存储一次
num_time_points = length(time_points);
time_values = time_points * dt;    % 实际时间值

% 寻找扭曲态函数
function initial_phase = generate_twisted_state(n, twist_number)
    % 生成扭曲态的初始相位
    positions = (0:n-1)'; % 振荡器的位置索引
    % 相位随位置线性增加，总增量为 2*pi*twist_number
    initial_phase = mod(2*pi*twist_number*positions/n, 2*pi) - pi;
end

% 初始化存储变量
twist_states = 0:3;  % 0，1，2，3扭曲态
num_states = length(twist_states);
p_order_history_all = zeros(num_states, num_time_points);
theta_history_all = zeros(n, num_time_points, num_states);

fprintf('模拟参数:\n');
fprintf('振荡器数量: %d\n', n);
fprintf('耦合范围: %d\n', r);
fprintf('成对耦合强度: %.1f\n', sigma);
fprintf('三元耦合强度: %.1f\n', sigma_delta);
fprintf('噪声强度: %.1f\n', D);
fprintf('总模拟时间: %d秒\n', T);

% 对每个扭曲态进行模拟
for state_idx = 1:num_states
    twist_number = twist_states(state_idx);
    fprintf('\n开始模拟 %d-扭曲态...\n', twist_number);
    
    % 生成对应扭曲态的初始相位
    initial_phase = generate_twisted_state(n, twist_number);
    
    % 初始化存储变量
    theta_history = zeros(n, num_time_points);
    p_order_history = zeros(1, num_time_points);
    
    % 初始化相位
    theta = initial_phase;
    time_point_idx = 1;
    
    % 运行模拟
    for t = 1:steps
        % 记录指定时间点的状态
        if ismember(t, time_points)
            theta_history(:, time_point_idx) = theta;
            
            % 计算局部序参量
            local_order = calculate_local_order(theta, n, r);
            
            % 计算P_order
            p_order_history(time_point_idx) = calculate_p_order(local_order);
            
            time_point_idx = time_point_idx + 1;
        end
        
        % 计算相互作用力
        dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
        
        % 更新相位
        noise = D * sqrt(dt) * randn(n, 1);
        theta = theta + dtheta * dt + noise;
        
        % 确保相位在[-π, π]范围内
        theta = mod(theta + pi, 2*pi) - pi;
    end
    
    % 存储当前扭曲态的结果
    p_order_history_all(state_idx, :) = p_order_history;
    theta_history_all(:, :, state_idx) = theta_history;
    
    fprintf('%d-扭曲态模拟完成\n', twist_number);
end

% 计算并显示总运行时间
total_time = toc;
fprintf('\n所有模拟完成，总运行时间: %.2f秒 (%.2f分钟)\n', total_time, total_time/60);

% 创建一行四列的子图
figure('Position', [50, 100, 1800, 400], 'Color', 'white');

% 创建颜色方案
colors = [
    0.9290, 0.6940, 0.1250; % 黄色 (0-twisted)
    0.4940, 0.1840, 0.5560; % 紫色 (1-twisted)
    0.3010, 0.7450, 0.9330; % 天蓝色 (2-twisted)
    0.8500, 0.3250, 0.0980; % 橙红色 (3-twisted)
];

% 创建图例句柄数组
legend_handles = zeros(num_states, 1);

% 设置统一的Y轴范围
y_lim = [0.8, 1.1];

% 对每个扭曲态创建一个子图
for i = 1:num_states
    twist_number = twist_states(i);
    
    % 创建子图
    ax = subplot(1, 4, i);
    
    % 获取当前扭曲态的数据
    p_order_data = p_order_history_all(i, :);
    
    % 绘制P_order曲线
    h = plot(time_values, p_order_data, '-o', 'Color', colors(i,:), 'LineWidth', 2.0, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'MarkerIndices', 1:50:length(time_values));
    
    % 存储句柄用于图例
    legend_handles(i) = h;
    
    % 设置子图标题（只包含扭曲态名称）
    title([num2str(twist_number), '-twisted state'], 'FontSize', 14);
    
    % 设置坐标轴标签
    xlabel('Time (s)', 'FontSize', 12);
    if i == 1
        ylabel('P_{order}', 'FontSize', 12, 'Interpreter', 'tex');
    end
    
    % 设置网格和边框
    grid on;
    set(gca, 'GridAlpha', 0.15, 'GridLineStyle', ':');
    box on;
    ylim(y_lim);  % 统一的Y轴范围
    
    % 设置轴和刻度的样式
    set(gca, 'TickDir', 'out');
    set(gca, 'LineWidth', 1.2);
    set(gca, 'FontSize', 12);
    
    % 添加最终P_order值
    final_p_order = p_order_data(end);
    text(max(time_values)*0.6, 0.6, {
        ['P_{final} = ', num2str(final_p_order, '%.3f')]
    }, 'FontSize', 12, 'BackgroundColor', [1 1 1 0.7]);
    
    % 创建更紧凑的布局，减少右侧空间
    pos = get(ax, 'Position');
    pos(3) = pos(3) * 0.95;  % 减少宽度，为图例腾出空间
    set(ax, 'Position', pos);
end

% 添加全局标题（只显示参数信息，不使用"标题"字样）
annotation('textbox', [0.35, 0.94, 0.3, 0.05], 'String', ...
    ['P_{order} Evolution for States (\sigma_{\Delta} = ', num2str(sigma_delta), ', D = ', num2str(D), ')'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

% 添加统一的图例到图表右侧
legend_labels = {'0-twisted', '1-twisted', '2-twisted', '3-twisted'};
h_legend = legend(legend_handles, legend_labels, 'Location', 'eastoutside', 'FontSize', 12);
set(h_legend, 'Position', [0.92, 0.3, 0.05, 0.4]);

% 调整整体图表大小，减少子图之间的间距
set(gcf, 'Units', 'pixels');
fig_pos = get(gcf, 'Position');
fig_pos(3) = 1650;  % 减小图表总宽度
set(gcf, 'Position', fig_pos);

% 保存结果
save('twisted_states_porder_evolution.mat', 'time_values', 'p_order_history_all', ...
    'theta_history_all', 'twist_states', 'sigma_delta', 'D', 'n', 'r', 'sigma');

% 辅助函数1: 计算相互作用
function dtheta = compute_interactions(theta, n, r, sigma, sigma_delta)
    dtheta = zeros(n, 1);
    
    % 计算所有振荡器的相位变化
    for i = 1:n
        % 计算成对相互作用项
        pair_interaction = 0;
        for j = i-r:i+r
            % 处理周期性边界条件
            j_idx = mod(j-1, n) + 1;
            if j_idx ~= i
                pair_interaction = pair_interaction + sin(theta(j_idx) - theta(i));
            end
        end
        pair_interaction = sigma * pair_interaction / (2*r);
        
        % 计算三元相互作用项
        triad_interaction = 0;
        count = 0;
        for j = i-r:i+r
            for k = i-r:i+r
                % 处理周期性边界条件
                j_idx = mod(j-1, n) + 1;
                k_idx = mod(k-1, n) + 1;
                % 确保三角形涉及三个不同的节点
                if j_idx ~= i && k_idx ~= i && j_idx ~= k_idx
                    triad_interaction = triad_interaction + sin(theta(j_idx) + theta(k_idx) - 2*theta(i));
                    count = count + 1;
                end
            end
        end
        if count > 0
            triad_interaction = sigma_delta * triad_interaction / (2*r*(2*r-1));
        end
        
        % 合并两种相互作用
        dtheta(i) = pair_interaction + triad_interaction;
    end
end

% 辅助函数2: 计算局部序参量
function local_order = calculate_local_order(theta, n, r)
    local_order = zeros(n, 1);
    
    for j = 1:n
        % 计算每个振荡器的局部序参量
        local_sum = 0;
        for k = j-r:j+r
            % 处理周期性边界条件
            k_idx = mod(k-1, n) + 1;
            local_sum = local_sum + exp(1i * theta(k_idx));
        end
        local_order(j) = abs(local_sum/(2*r+1));
    end
end

% 辅助函数3: 计算P_order
function p_order = calculate_p_order(local_order)
    % 计算有序振荡器比例 - 使用阈值0.8
    disordered_count = sum(local_order < 0.8);
    p_order = 1 - disordered_count/length(local_order);
end