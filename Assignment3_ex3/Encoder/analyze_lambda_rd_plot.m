function analyze_lambda_rd_plot()
    % Basic setup
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
    
    % Fixed QP
    QP = 1;
    
    % More fine-grained lambda values
    lambda_values = [0.2:0.1:1.2];  % Smaller steps for smoother curve
    
    % Initialize results arrays
    total_bits = zeros(1, length(lambda_values));
    psnrs = zeros(1, length(lambda_values));
    
    % Pre-process video
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);
    
    % Test each lambda value
    for i = 1:length(lambda_values)
        lambda = lambda_values(i);
        fprintf('Testing lambda = %.2f\n', lambda);
        
        % Run encoder with VBS enabled
        encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, ...
               blockSize, searchRange, blockSize, QP, I_Period, 1, ...
               lambda, true, false, false);
        
        % Run decoder and get results
        [total_bytes, ~] = decoder(decodedFile);
        psnr_val = calculatePSNR(decodedFile, outputFile, width, height, numFrames);
        
        % Store results
        total_bits(i) = total_bytes * 8;
        psnrs(i) = psnr_val;
    end
    
    % Create R-D plot
    figure;
    plot(total_bits, psnrs, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    grid on;
    xlabel('totalBits');
    ylabel('PSNR');
    title('R-D plot when QP=1 varying Lambda');
    
    % Format axes
    ax = gca;
    ax.XAxis.Exponent = 6;  % Show x-axis in millions
    
    % Add data points labels if desired
    for i = 1:length(lambda_values)
        text(total_bits(i), psnrs(i), sprintf(' λ=%.1f', lambda_values(i)), ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 8);
    end
    
    % Save the plot
    saveas(gcf, '../Outputs/rd_plot_qp1.png');
    
    % Print results
    fprintf('\nResults:\n');
    fprintf('Lambda\tPSNR\tBits\n');
    fprintf('------------------------\n');
    for i = 1:length(lambda_values)
        fprintf('%.2f\t%.2f\t%d\n', lambda_values(i), psnrs(i), total_bits(i));
    end
end