clear; close all; clc;

addpath('../Utils');    % For utils functions
addpath('../Outputs');  % For output files handling
addpath('../Decoder');  % For Decoder functions
addpath('../Encoder');  % For Encoder functions

% The input file: 21 frames total = 7 Foreman CIF + 7 Akiyo CIF + 7 Foreman CIF frames
filename = '../CIF.yuv';
outputFile = '../Outputs/Y_only_foreman.yuv';
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv';
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';

% CIF dimensions
width = 352;
height = 288;

numFrames = 21;
blockSize = 8;
searchRange = 16;
I_Period = 21; % As per exercise requirement

% Bitcount estimation arrays for CIF
i_bitCountPerRow = [25051, 25360, 18853, 12997, 8962, 5237, 3086, 1714, 737, 439,389,391];
p_bitCountPerRow = [19704, 19646,14576,9974,6222,3538,2058,1293,734,475,378,354];

fps = 30;

% QPs and Bitrates
QPs = [3, 6, 9];
bitrates = [7983360, 2737152, 410572]; % in bits per second
%bitrates = [2737152]; % in bits per second
% bitrates(1) ~7.98 Mbps, bitrates(2) ~2.737 Mbps, bitrates(3) ~0.411 Mbps

% Build configuration array
configs = [];

% RCflag=0, vary QP
rc = 0;
for q = 1:length(QPs)
    cfg.name = sprintf('RC_%d_QP_%d', rc, QPs(q));
    cfg.RCflag = rc;
    cfg.QP = QPs(q);
    cfg.targetBR = 0; % Not used for RCflag=0
    cfg.i_period = I_Period;
    configs = [configs, cfg];
end

% RCflag=1,2 vary bitrate
for rc = 1:2
    for b = 1:length(bitrates)
        cfg.name = sprintf('RC_%d_BR_%0.1fk', rc, bitrates(b)/1000);
        cfg.RCflag = rc;
        cfg.QP = 0; % QP not fixed for rate control mode
        cfg.targetBR = bitrates(b);
        cfg.i_period = I_Period;
        configs = [configs, cfg];
    end
end

% Create output directory if it doesn't exist
if ~exist('../Outputs', 'dir')
    mkdir('../Outputs');
end

% Initialize results storage
results = struct();
for i = 1:length(configs)
    results(i).name = configs(i).name;
    results(i).bitrate = 0;
    results(i).psnr = 0;
    results(i).psnrValues = [];  % <-- Store per-frame PSNR
    results(i).time = 0;
    results(i).RCflag = configs(i).RCflag;
    results(i).QP = configs(i).QP;
    results(i).targetBR = configs(i).targetBR;
end

% Encode and decode for each configuration
for c = 1:length(configs)
    fprintf('Testing configuration: %s\n', configs(c).name);
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

    RCflag = configs(c).RCflag;
    QP = configs(c).QP;
    targetBR = configs(c).targetBR;
    i_period = configs(c).i_period;
    per_frame_budget = 0;
    if (RCflag > 1)
        per_block_row_budget = get_per_row_budget(targetBR,fps,paddedWidth,paddedHeight,blockSize);
        QP = findCorrectQP(per_block_row_budget,p_bitCountPerRow);
        QP = QP - 1;
    end

    if RCflag == 0
        % Constant QP mode
        per_block_row_budget = 0; 
        RCflag_val = false;
    else
        [per_block_row_budget,per_frame_budget] = get_per_row_budget(targetBR, fps, width, height, blockSize);
        
        RCflag_val = RCflag;
    end

    lambda = get_lambda_for_qp(QP);

    % Time encoding + decoding
    tstart = tic;
    encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, ...
            blockSize, searchRange, blockSize, QP, i_period, 1, ...
            lambda, true, true, true, RCflag_val, per_block_row_budget, ...
            p_bitCountPerRow, i_bitCountPerRow,per_frame_budget);

    [total_bytes, total_bytes_per_frame] = decoder(decodedFile);
    time_elapsed = toc(tstart);

    [psnr_val, psnrValues] = calculatePSNR(decodedFile, outputFile, width, height, numFrames);

    % Store results
    results(c).bitrate = total_bytes;  % bits/sec
    results(c).psnr = psnr_val;
    results(c).psnrValues = psnrValues;
    results(c).time = time_elapsed;
end

% Display results
for i = 1:length(configs)
    fprintf('Configuration: %s\n', configs(i).name);
    fprintf('  Rate Control: %d\n', configs(i).RCflag);
    fprintf('  Target Bitrate: %.2f kbps\n', results(i).targetBR / 1000);
    fprintf('  I-Period: %d\n', configs(i).i_period);
    fprintf('  PSNR: %.2f dB\n', results(i).psnr);
    fprintf('  Achieved Bitrate: %.2f bits/s\n', results(i).bitrate);
    fprintf('\n');
end

% Separate results by RCflag
rc0_indices = find([results.RCflag] == 0);
rc1_indices = find([results.RCflag] == 1);
rc2_indices = find([results.RCflag] == 2);

% Sort them according to QP or bitrate
[~,idx_rc0] = sort([results(rc0_indices).QP]);
rc0_results = results(rc0_indices(idx_rc0));

[~,idx_rc1] = sort([results(rc1_indices).targetBR],'descend');
rc1_results = results(rc1_indices(idx_rc1));

[~,idx_rc2] = sort([results(rc2_indices).targetBR],'descend');
rc2_results = results(rc2_indices(idx_rc2));

% Plot R-D curves
figure;
hold on;
plot([rc0_results.bitrate], [rc0_results.psnr], '-o', 'DisplayName','RCflag=0');
plot([rc1_results.bitrate], [rc1_results.psnr], '-x', 'DisplayName','RCflag=1');
plot([rc2_results.bitrate], [rc2_results.psnr], '-s', 'DisplayName','RCflag=2');
xlabel('Bitrate (bps)');
ylabel('PSNR (dB)');
legend('Location','best');
title('R-D Curves for Different RCflag Modes');
grid on;

% Create a table of encoding times
all_names = {results.name}';
all_times = [results.time]';
TimeTable = table(all_names, all_times, 'VariableNames', {'Experiment','EncodingTime_s'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot per-frame PSNR for RCflag=1,2 at ~2.4 Mbps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We know bitrates(2) = 2737152 ~ 2.4 Mbps
target_bitrate_2_4mbps = 2737152;

% Find results for RCflag=1,2 and targetBR ~2.4Mbps
rc1_2_4mbps_idx = find([results.RCflag] == 1 & [results.targetBR] == target_bitrate_2_4mbps);
rc2_2_4mbps_idx = find([results.RCflag] == 2 & [results.targetBR] == target_bitrate_2_4mbps);

% NOTE: Original instructions mentioned RCflag=1,2,3 for 2 Mbps, but we only have RCflag=1,2
% If you had RCflag=3 experiments as well, just find them similarly:
rc3_2_4mbps_idx = find([results.RCflag] == 3 & [results.targetBR] == target_bitrate_2_4mbps);

% Ensure these indices exist. If you didn't run RCflag=3, just comment that section out.
if isempty(rc3_2_4mbps_idx)
    warning('No RCflag=3 configuration found for 2.4 Mbps. The third curve will not be plotted.');
end

figure; hold on;
if ~isempty(rc1_2_4mbps_idx)
    plot(results(rc1_2_4mbps_idx).psnrValues, '-o', 'DisplayName', 'RCflag=1');
end
if ~isempty(rc2_2_4mbps_idx)
    plot(results(rc2_2_4mbps_idx).psnrValues, '-x', 'DisplayName', 'RCflag=2');
end
% if ~isempty(rc3_2_4mbps_idx)
%     plot(results(rc3_2_4mbps_idx).psnrValues, '-s', 'DisplayName', 'RCflag=3');
% end

xlabel('Frame Number');
ylabel('PSNR (dB)');
title('Per-frame PSNR at ~2.4 Mbps for RCflag=1,2,3');
legend('Location','best');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ Lambda function ------
function lambda = get_lambda_for_qp(QP)
    if QP == 0
        lambda = 0.01;
    elseif QP == 1
        lambda = 0.02;
    elseif QP == 2
        lambda = 0.1667;
    elseif QP == 3
        lambda = 0.1667;
    elseif QP == 4
        lambda = 0.1667;
    elseif QP == 5
        lambda = 0.1852;
    elseif QP == 6
        lambda = 0.2037;
    elseif QP == 7
        lambda = 0.222;
    elseif QP == 8
        lambda = 0.3148;
    elseif QP == 9
        lambda = 0.4074;
    elseif QP == 10
        lambda = 0.5;
    elseif QP == 11
        lambda = 0.55;
    else
        lambda = NaN;
    end
end
