function analyze_reference_frames()
    % Parameters
      filename = '../CIF.yuv';  % YUV file to read
     % filename = '../QCIF.yuv';  % YUV file to read
    outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
    paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
    referenceFile = '../Outputs/referenceFrames.yuv';
    decodedFile = '../Outputs/decoded_Y_foreman.yuv';
    width = 352;                     % Frame width
    height = 288;                    % Frame height
    % width = 176;                     % Frame width
    % height = 144;                    % Frame height
    numFrames = 21;                 % Number of frames tfalseo process
    searchRange = 16;                 % Search range r = 1,4, and 8
    QP = 7;
    j = 4;
    VBSEnable = true;
    FMEEnable = true;
    FastME = true;
    RCflag = 0;
    
    if(VBSEnable == true)
            j = j-1;
    end
    
    blockSize = 2^j;                   % Block size for motion estimation
    dct_blockSize = 2^j;
    I_Period = 21; 
    nRefFrames = 1;                 % Can take value from 1 to 4
    
      targetBR = 2737152;
      % targetBR = 1094860;
    
    
    fps = 30;
    %bitCountPerRow = [2112, 1520, 1256, 1156, 1076, 1056, 65, 32, 15, 7];
  
% For CIF
i_bitCountPerRow = [25051, 25360, 18853, 12997, 8962, 5237, 3086, 1714, 737, 439,389,391];
p_bitCountPerRow = [19704, 19646,14576,9974,6222,3538,2058,1293,734,475,378,354];
% For QCIF
% i_bitCountPerRow = [13599,13781,10479,7783,5642,3870,2479,1326,489,216,186,181];
% p_bitCountPerRow = [9608, 9624,7059,4936,3315,2121,1318,824,433,253,189,183];

    
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
            lambda = 0.2037
        elseif QP == 7
            lambda = 0.222;
        elseif QP == 8
            lambda = 0.3148;
        elseif QP == 9
            lambda = 0.4074
        elseif QP == 10
            lambda = 0.5;
        elseif QP == 11
            lambda = 0.55;
        else
            % Handle other cases if needed
            lambda = NaN;  % Example: set to NaN if QP is not 1, 4, 7, or 10
       end
    end
    lambda = get_lambda_for_qp(QP);

    % Initialize storage for results
    max_ref_frames = 1;
    psnrs = zeros(max_ref_frames, numFrames);
    bits_per_frame = zeros(max_ref_frames, numFrames);


    % Create output directory if it doesn't exist
    if ~exist('../Outputs', 'dir')
        mkdir('../Outputs');
    end

    % Pre-process video
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

    per_block_row_budget = 0;
    if RCflag == true
        per_block_row_budget = get_per_row_budget(targetBR,fps,paddedWidth,paddedHeight,blockSize);
    end
    % Test each reference frame configuration
    for nRefFrames = 1:max_ref_frames
        fprintf('Testing with %d reference frame(s)\n', nRefFrames);

        % Run encoder
        encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME,RCflag,per_block_row_budget,p_bitCountPerRow,i_bitCountPerRow);

        % Run decoder and collect frame-by-frame statistics
        [total_byte, bytes_list] = decoder(decodedFile);
        bits_per_frame(nRefFrames, :) = bytes_list ; 

        % Calculate PSNR for each frame
        fidOriginal = fopen(outputFile, 'r');
        fidDecoded = fopen(decodedFile, 'r');

        for frame = 1:numFrames
            % Read original frame
            originalFrame = fread(fidOriginal, [width, height], 'uint8')';
            % Read decoded frame
            decodedFrame = fread(fidDecoded, [width, height], 'uint8')';

            % Calculate MSE and PSNR
            mse = mean((double(originalFrame(:)) - double(decodedFrame(:))).^2);
            if mse > 0
                psnrs(nRefFrames, frame) = 10 * log10(255^2 / mse);
            else
                psnrs(nRefFrames, frame) = Inf;
            end
        end

        fclose(fidOriginal);
        fclose(fidDecoded);
    end

    % Create plots
    figure('Position', [100, 100, 1200, 500]);

    % Plot PSNR per frame
    subplot(1, 2, 1);
    markers = {'o-', 's-', 'd-', '^-'};
    for ref = 1:max_ref_frames
        plot(1:numFrames, psnrs(ref, :), markers{ref}, 'LineWidth', 2, ...
             'DisplayName', sprintf('%d Ref Frame(s)', ref));
        hold on;
    end
    grid on;
    xlabel('Frame Number');
    ylabel('PSNR (dB)');
    title(sprintf('PSNR per Frame for CIF I Period = %d',I_Period));
    % legend('Location', 'best');
    legend off;
    hold off;

    % Plot bits per frame
    subplot(1, 2, 2);
    for ref = 1:max_ref_frames
        plot(1:numFrames, bits_per_frame(ref, :), markers{ref}, 'LineWidth', 2, ...
             'DisplayName', sprintf('%d Ref Frame(s)', ref));
        hold on;
    end
    grid on;
    xlabel('Frame Number');
    ylabel('Frame Size (bits)');
    title(sprintf('Bits per Frame for CIF I Period = %d', I_Period));
    % legend('Location', 'best');
    legend off;
    hold off;

    % Save the plots
    saveas(gcf, '../Outputs/reference_frames_analysis.png');

    % Print summary statistics
    fprintf('\nSummary Statistics:\n');
    fprintf('==================\n');
    for ref = 1:max_ref_frames
        fprintf('\nReference Frames: %d\n', ref);
        fprintf('Average PSNR: %.2f dB\n', mean(psnrs(ref, :)));
        fprintf('Average bits per frame: %.2f\n', mean(bits_per_frame(ref, :)));
        fprintf('Total bits: %.2f\n', sum(bits_per_frame(ref, :)));
    end

    % Save numerical results
    save('../Outputs/reference_frames_results.mat', 'psnrs', 'bits_per_frame');
end