function generate_rd_analysis()
    % Fixed parameters
    filename = '../foreman_cif-1.yuv';
    outputFile = '../Outputs/Y_only_foreman.yuv';
    paddedOutputFile = '../Outputs/padded_Y_foreman.yuv';
    referenceFile = '../Outputs/referenceFrames.yuv';
    decodedFile = '../Outputs/decoded_Y_foreman.yuv';
    width = 352;
    height = 288;
    numFrames = 10;
    blockSize = 16;
    searchRange = 4;
    I_Period = 8;
    QPs = [1, 4, 7, 10];
    lambda = 1;

    % Define configurations as a struct array
    configs = struct([]);
    
    % Base encoder
    configs(1).vbs = false;
    configs(1).fme = false;
    configs(1).fastme = false;
    configs(1).nref = 1;
    configs(1).name = 'Base';
    
    % VBS only
    configs(2).vbs = true;
    configs(2).fme = false;
    configs(2).fastme = false;
    configs(2).nref = 1;
    configs(2).name = 'VBS only';
    
    % FME only
    configs(3).vbs = false;
    configs(3).fme = true;
    configs(3).fastme = false;
    configs(3).nref = 1;
    configs(3).name = 'FME only';
    
    % FastME only
    configs(4).vbs = false;
    configs(4).fme = false;
    configs(4).fastme = true;
    configs(4).nref = 1;
    configs(4).name = 'FastME only';
    
    % Multiple reference frames only
    configs(5).vbs = false;
    configs(5).fme = false;
    configs(5).fastme = false;
    configs(5).nref = 4;
    configs(5).name = 'MultiRef only';
    
    % All features
    configs(6).vbs = true;
    configs(6).fme = true;
    configs(6).fastme = true;
    configs(6).nref = 4;
    configs(6).name = 'All Features';

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
            
            % Run encoder with current configuration
            encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, ...
                   blockSize, searchRange, blockSize, QP, I_Period, configs(c).nref, ...
                   lambda, configs(c).vbs, configs(c).fme, configs(c).fastme);
            
            % Run decoder
            [total_bytes, bytes_list] = decoder(decodedFile);
            
            % Store execution time
            results(c).times(qp_idx) = toc;

            % Calculate PSNR
            psnr_val = calculatePSNR(decodedFile, outputFile, width, height, numFrames);
            
            % Store results (convert bytes to bits)
            results(c).bitrates(qp_idx) = total_bytes * 8; % Convert to bits per frame
            results(c).psnrs(qp_idx) = psnr_val;
        end
    end

    % Create figure for plots
    fig = figure('Position', [100, 100, 1200, 500]);

    % RD curve subplot
    subplot(1, 2, 1);
    colors = lines(length(configs));
    markers = {'o', 's', 'd', '^', 'v', 'p'};
    
    for i = 1:length(results)
        plot(results(i).bitrates, results(i).psnrs, ...
             ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', results(i).name, 'MarkerSize', 8);
        hold on;
    end
    
    grid on;
    xlabel('Total Size in Bits');
    ylabel('PSNR (dB)');
    title('Rate-Distortion Curve');
    legend('Location', 'southeast');
    ax = gca;
    ax.XAxis.Exponent = 6;  % Force x10^6 scaling
    xlim([0 4.5e6]);        % Match the reference graph range
    ylim([15 50]);          % Match the reference graph range
    hold off;

    % Execution time subplot
    subplot(1, 2, 2);
    for i = 1:length(results)
        plot(QPs, results(i).times, ...
             ['-' markers{i}], 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', results(i).name, 'MarkerSize', 8);
        hold on;
    end
    
    grid on;
    xlabel('QP');
    ylabel('Execution Time (seconds)');
    title('Execution Time vs QP');
    legend('Location', 'northwest');
    hold off;

    % Save the figure
    saveas(fig, '../Outputs/rd_analysis.png');
    
    % Save numerical results
    save('../Outputs/rd_results.mat', 'results', 'QPs');
    
    % Print summary statistics
    print_summary_statistics(results, QPs);
end

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