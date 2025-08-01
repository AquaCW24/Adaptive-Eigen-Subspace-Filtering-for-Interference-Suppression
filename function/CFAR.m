function [INDEX, THRESHOLD] = CFAR(INPUT, TRAINING_CELL, GUARD_CELL, PFA, MODE, PLOT_BOOL)
% CFAR implementation for target detection
% Mode: "CA", "GOCA", "SOCA"
% Input:
% - INPUT: MxN "power" Matrix
% - TRAINING_CELL: training cell, number of single side
% - GUARD_CELL: guard cell, number of single side
% - PFA: Probability of false alarm
% - MODE: The metric of noise estimation ("CA", "GOCA", "SOCA")
% - PLOT_BOOL: 
% 
% Outputs: 
% INDEX: Indices of detected targets
% THRESHOLD: Threshold for detecting the appearance of a target

data = mean(INPUT,1);
Std = std(10*log10(INPUT),0,1);
threshold_std = mean(Std) - 2*std(Std);

N = length(data); % Signal length
alpha = 2 * TRAINING_CELL * (PFA^(-1/(2*TRAINING_CELL)) - 1);  % Threshold factor

THRESHOLD = zeros(1, N);

% Process each cell with sliding window
for i = 1:N
    % Define window boundaries
    left_start = max(1, i - TRAINING_CELL - GUARD_CELL);
    left_end = max(1, i - GUARD_CELL - 1);
    right_start = min(N, i + GUARD_CELL + 1);
    right_end = min(N, i + TRAINING_CELL + GUARD_CELL);

    % Extract training cells
    left_window = data(left_start:left_end);
    right_window = data(right_start:right_end);

    if i <= (GUARD_CELL + 1)
        left_window = [];
    elseif i >= (N - GUARD_CELL - 1)
        right_window = [];
    end
    
    % Calculate average cell noise power
    if MODE == "CA" % Cell Average
        noise_power = (sum(left_window) + sum(right_window)) / ...
                    (length(left_window) + length(right_window));
    elseif MODE == "GOCA" % Greatest-of CFAR
        noise_power = max(mean(left_window), mean(right_window));
    elseif MODE == "SOCA" % Smallest Of CFAR
        noise_power = min(mean(left_window), mean(right_window));
    elseif MODE == "OS" % Order Statistic
        noise_power = sort([left_window, right_window]);
        noise_power = noise_power(floor(0.75*length([left_window, right_window])));
    end

    % Calculate threshold
    THRESHOLD(i) = alpha * noise_power;

end

% Obtain target index
INDEX = find(data > THRESHOLD | Std < threshold_std);

if PLOT_BOOL == "true"
    % Plotting target
    figure('Units','normalized','Position',[0 .1 .9 .4]);
    plot(10*log10(data), 'k', 'LineWidth', 1);
    hold on;
    plot(10*log10(THRESHOLD), 'r', 'LineWidth', 1);
    plot(INDEX, 10*log10(data(INDEX)), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    grid on;box on;
    title({MODE+'-CFAR Detection';sprintf('Training Cell: %d, Guard Cell: %d, P_{FA}: %.3f, P_{FA(actual)}:%.3f',TRAINING_CELL*2,GUARD_CELL,PFA,length(INDEX)/N)});
    xlabel('Samples');
    ylabel('Power (dB)');
    legend('Signal', 'Threshold', 'Detected Targets');
    xlim([1 N])
end

end