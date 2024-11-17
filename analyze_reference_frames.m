function analyze_reference_frames()
    % Fixed parameters
    filename = '../synthetic.yuv';
    outputFile = '../Outputs/Y_only_synthetic.yuv';
    paddedOutputFile = '../Outputs/padded_Y_synthetic.yuv';
    referenceFile = '../Outputs/referenceFrames.yuv';
    decodedFile = '../Outputs/decoded_Y_synthetic.yuv';
    width = 352;
    height = 288;
    numFrames = 10;
    blockSize = 16;
    searchRange = 4;
    I_Period = 8;
    QP = 4;
    lambda = 1;

    % Initialize storage for results
    max_ref_frames = 4;
    psnrs = zeros(max_ref_frames, numFrames);
    bits_per_frame = zeros(max_ref_frames, numFrames);

    % Base configuration (only vary nRefFrames)
    VBSEnable = false;
    FMEEnable = false;
    FastME = false;

    % Create output directory if it doesn't exist
    if ~exist('../Outputs', 'dir')
        mkdir('../Outputs');
    end

    % Pre-process video
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

    % Test each reference frame configuration
    for nRefFrames = 1:max_ref_frames
        fprintf('Testing with %d reference frame(s)\n', nRefFrames);

        % Run encoder
        encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, ...
               blockSize, searchRange, blockSize, QP, I_Period, nRefFrames, ...
               lambda, VBSEnable, FMEEnable, FastME);

        % Run decoder and collect frame-by-frame statistics
        [~, bytes_list] = decoder(decodedFile);
        bits_per_frame(nRefFrames, :) = bytes_list * 8; % Convert to bits

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
    title('PSNR per Frame for Different Numbers of Reference Frames');
    legend('Location', 'best');
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
    title('Bits per Frame for Different Numbers of Reference Frames');
    legend('Location', 'best');
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