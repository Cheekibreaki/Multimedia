function analyze_lambda_qp_relation()
    % Basic setup (use your existing parameters)
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
    
    % Choose two representative QP values
    QPs = 1;  % q = 2 candidate QP values
    
    % Lambda values to test for each QP
    lambda_values = [0.2:0.1:1.2];
    
    % Initialize results storage
    results = zeros(length(QPs), length(lambda_values), 2); % Store PSNR and Bitrate
    
    % Test each combination
    for qp_idx = 1:length(QPs)
        for lambda_idx = 1:length(lambda_values)
            QP = QPs(qp_idx);
            lambda = lambda_values(lambda_idx);
            
            % Run encoder with VBS enabled
            encoder(referenceFile, paddedOutputFile, numFrames, width, height, ...
                   blockSize, searchRange, blockSize, QP, I_Period, 1, ...
                   lambda, true, false, false);
            
            % Run decoder and get results
            [total_bytes, ~] = decoder(decodedFile);
            psnr_val = calculatePSNR(decodedFile, outputFile, width, height, numFrames);
            
            % Store results
            results(qp_idx, lambda_idx, 1) = psnr_val;
            results(qp_idx, lambda_idx, 2) = total_bytes * 8;  % Convert to bits
            
            fprintf('QP=%d, lambda=%.2f: PSNR=%.2f, Bits=%d\n', ...
                    QP, lambda, psnr_val, total_bytes * 8);
        end
    end
    
    % Plot results
    figure;
    for qp_idx = 1:length(QPs)
        subplot(1, 2, qp_idx);
        plot(lambda_values, results(qp_idx, :, 1), '-o');
        xlabel('Lambda');
        ylabel('PSNR (dB)');
        title(sprintf('PSNR vs Lambda for QP=%d', QPs(qp_idx)));
        grid on;
    end
    
    % Find best lambda for each QP
    best_lambdas = zeros(1, length(QPs));
    for qp_idx = 1:length(QPs)
        % Find lambda that gives best RD performance
        rd_scores = results(qp_idx, :, 1) ./ log(results(qp_idx, :, 2));
        [~, best_idx] = max(rd_scores);
        best_lambdas(qp_idx) = lambda_values(best_idx);
    end
    
    % Create lambda-QP relation
    p = polyfit(QPs, best_lambdas, 1);  % Linear fit
    fprintf('\nSuggested lambda-QP relation:\n');
    fprintf('lambda = %.3f * QP + %.3f\n', p(1), p(2));
    
    % Plot lambda-QP relation
    figure;
    plot(QPs, best_lambdas, 'o', 'DisplayName', 'Best lambdas');
    hold on;
    qp_range = [min(QPs)-1, max(QPs)+1];
    plot(qp_range, polyval(p, qp_range), '-', 'DisplayName', 'Fitted relation');
    xlabel('QP');
    ylabel('Best Lambda');
    title('Lambda-QP Relation');
    legend;
    grid on;
end