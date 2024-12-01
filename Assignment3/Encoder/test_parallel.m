addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../CIF.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 20;                 % Number of frames to process
searchRange = 1;                 % Search range r = 1,4, and 8
j = 4;
VBSEnable = false;
FMEEnable = false;
FastME = false;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 8; 
nRefFrames = 4;                 % Can take value from 1 to 4

QPs = [1, 2, 4, 7, 10];
lambdas = [0,0.166667,0.166667,0.22222,0.5];

coreCount = 2;
% Start parallel pool
parpool(coreCount);

    % Define configurations as a struct array
    configs = struct([]);
    
    % Mode 0, original encoder, no parallism
    configs(1).mode = 0;
    configs(1).name = 'Mode 0';

    % Mode 1, extreme block-level parallelism
    configs(2).mode = 1;
    configs(2).name = 'Mode 1';

    % Mode 2
    configs(3).mode = 2;
    configs(3).name = 'Mode 2';

    % Mode 3 
    configs(4).mode = 3;
    configs(4).name = 'Mode 3';

    % Initialize results storage
    results = struct();
    for i = 1:length(configs)
        results(i).name = configs(i).name;
        results(i).bitrates = zeros(1, length(QPs));
        results(i).psnrs = zeros(1, length(QPs));
        results(i).times = zeros(1, length(QPs));
    end

        % Create output directory if it doesn't exist
    if ~exist('../Outputs', 'dir')
        mkdir('../Outputs');
    end

    colors = lines(length(configs));

    % Pre-process video
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

    % Run each configuration
    for c = 1:length(configs)
        fprintf('Testing configuration: %s\n', configs(c).name);
        
        for qp_idx = 1:length(QPs)
            QP = QPs(qp_idx);
            fprintf('  QP = %d\n', QP);

            % Time the encoding and decoding process
            tic;
                
            if configs(c).mode == 1
                %for mode 1 
                encoder_mode1(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames, lambdas(qp_idx), VBSEnable, FMEEnable, FastME, configs(c).mode);
            elseif configs(c).mode == 3
                encoder_mode3(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames, lambdas(qp_idx), VBSEnable, FMEEnable, FastME, configs(c).mode);
            else    
                %for mode 0 and mode 2
                encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames, lambdas(qp_idx), VBSEnable, FMEEnable, FastME,configs(c).mode);
            end
            % Run decoder
            [total_bytes, bytes_list] = decoder(decodedFile,configs(c).mode);
            
            % Store execution time
            results(c).times(qp_idx) = toc;

            % Calculate PSNR
            psnr_val = calculatePSNR(decodedFile, outputFile, width, height, numFrames);
            
            % Store results (convert bytes to bits)
            results(c).bitrates(qp_idx) = total_bytes * 8; % Convert to bits per frame
            results(c).psnrs(qp_idx) = psnr_val;
        end
    end
    % Delete the parallel pool
    delete(gcp('nocreate'));

    % Create figure for plots
    fig = figure('Position', [100, 100, 1200, 500]);

    % RD curve subplot
    subplot(2, 1, 1);
    
    markers = {'o', 's', 'd', '^', 'v', 'p', 'h', 'x'};
    
    for i = 1:length(results)
        plot(results(i).bitrates, results(i).psnrs, ...
             ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', results(i).name, 'MarkerSize', 8);
        hold on;
    end
    
    grid on;
    xlabel('Total Size in Bits');
    ylabel('PSNR (dB)');
    title('Rate-Distortion Curve for Different Parallel Modes');
    legend('Location', 'southeast');
    ax = gca;
    ax.XAxis.Exponent = 6;  % Force x10^6 scaling
    xlim([0 4.5e6]);        % Match the reference graph range
    ylim([15 50]);          % Match the reference graph range
    hold off;

    % Execution time subplot
    subplot(2, 1, 2);
    for i = 1:length(results)
        plot(QPs, results(i).times, ...
             ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', results(i).name, 'MarkerSize', 8);
        hold on;
    end
    
    grid on;
    xlabel('QP');
    ylabel('Execution Time (seconds)');
    title('Encoding Time for Different Parallel Modes');
    legend('Location', 'northwest');
    hold off;

    % Save the figure
    saveas(fig, '../Outputs/mode_compare.png');
    
    % Save numerical results
    save('../Outputs/rd_results.mat', 'results', 'QPs');
    
    % Print summary statistics
    print_summary_statistics(results, QPs);




function print_summary_statistics(results, QPs)
    fprintf('\nSummary Statistics:\n');
    fprintf('==================\n\n');
    
    % Calculate average improvements relative to base encoder
    base_idx = 1;  % Index of base encoder results
    
    for i = 2:length(results)
        % Calculate average PSNR improvement
        psnr_improvement = mean(results(i).psnrs - results(base_idx).psnrs);
        
        % Calculate average bitrate savings (%)
        bitrate_saving = mean((results(base_idx).bitrates - results(i).bitrates) ./ results(base_idx).bitrates * 100);
        
        % Calculate average time difference (%)
        time_diff = mean((results(i).times - results(base_idx).times) ./ results(base_idx).times * 100);
        
        fprintf('%s vs Base:\n', results(i).name);
        fprintf('  Average PSNR improvement: %.2f dB\n', psnr_improvement);
        fprintf('  Average bitrate saving: %.2f%%\n', bitrate_saving);
        fprintf('  Average time difference: %+.2f%%\n\n', time_diff);
    end
end



