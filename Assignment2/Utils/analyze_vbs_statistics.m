function analyze_vbs_statistics(filename, outputFile, paddedOutputFile, referenceFile, decodedFile, ...
                                 width, height, numFrames, blockSize, searchRange, ...
                                 I_Period, lambda, VBSEnable, FMEEnable, FastME, ...
                                 nRefFrames, QPs)

    % Initialize storage for VBS statistics
    split_percentages = zeros(1, length(QPs));
    bits_per_frame = zeros(1, length(QPs));

    % Configuration: Only VBS enabled
    VBSEnable = true;
    FMEEnable = false;
    FastME = false;
    nRefFrames = 1;

    % Create output directory if it doesn't exist
    if ~exist('../Outputs', 'dir')
        mkdir('../Outputs');
    end

    % Pre-process video
    fprintf('Pre-processing video...\n');
    dumpYComponentsToFile(filename, width, height, numFrames, outputFile);
    [paddedWidth, paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

    % Test each QP value
    for qp_idx = 1:length(QPs)
        QP = QPs(qp_idx);
        fprintf('Processing QP = %d\n', QP);

        % Run encoder
        encoder(referenceFile, paddedOutputFile, numFrames, paddedWidth, paddedHeight, ...
               blockSize, searchRange, blockSize, QP, I_Period, nRefFrames, ...
               lambda, VBSEnable, FMEEnable, FastME);

        % Run decoder to get total bytes
        [total_bytes, ~] = decoder(decodedFile);
        bits_per_frame(qp_idx) = total_bytes * 8;  % Convert bytes to bits

        % Calculate percentage of split blocks
        total_split_blocks = 0;
        total_blocks = 0;

        % Load and analyze VBS decisions from saved files
        for frameIdx = 1:numFrames
            % Skip I-frames when counting splits
            if mod(frameIdx-1, I_Period) ~= 0
                MDiffFile = sprintf('../Outputs/MDiff_frame_%d.mat', frameIdx);
                if ~exist(MDiffFile, 'file')
                    warning('Motion vector file not found: %s', MDiffFile);
                    continue;
                end

                % Load motion vector data
                load(MDiffFile, 'encodedMDiff');
                
                % Check if encodedMDiff exists and is not empty
                if ~exist('encodedMDiff', 'var') || isempty(encodedMDiff)
                    warning('Empty motion vector data in frame %d', frameIdx);
                    continue;
                end

                % Get binary decoded data
                motionVectorRevEGC = exp_golomb_decode(encodedMDiff(2:end));
                
                % Calculate expected number of 2x2 blocks
                num_2x2_blocks = (paddedHeight/blockSize/2) * (paddedWidth/blockSize/2);
                fprintf('Frame %d: Expected blocks = %d\n', frameIdx, num_2x2_blocks);

                % Process each 2x2 block
                idx = 1;
                block_count = 0;
                while idx <= length(motionVectorRevEGC) && block_count < num_2x2_blocks
                    if idx > length(motionVectorRevEGC)
                        break;
                    end

                    % Get the split flag (prefix)
                    prefix = motionVectorRevEGC(idx);
                    idx = idx + 1;
                    block_count = block_count + 1;

                    if prefix == 1  % Split block
                        total_split_blocks = total_split_blocks + 1;
                        idx = idx + 12;  % Skip split block motion vectors
                    else  % Non-split block
                        idx = idx + 3;   % Skip non-split block motion vector
                    end
                    total_blocks = total_blocks + 1;

                    fprintf('Frame %d, Block %d: prefix=%d\n', frameIdx, block_count, prefix);
                end
            end
        end

        % Calculate percentage of split blocks
        if total_blocks > 0
            split_percentages(qp_idx) = (total_split_blocks / total_blocks) * 100;
        else
            split_percentages(qp_idx) = 0;
            warning('No blocks processed for QP = %d', QP);
        end
        
        fprintf('QP=%d: Total blocks=%d, Split blocks=%d, Percentage=%.2f%%\n', ...
                QP, total_blocks, total_split_blocks, split_percentages(qp_idx));
    end

    % Create figures for analysis
    figure('Position', [100, 100, 1200, 500]);

    % Plot percentage of splits vs QP
    subplot(1, 2, 1);
    plot(QPs, split_percentages, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('QP Value');
    ylabel('Percentage of Split Blocks (%)');
    title('Block Split Decisions vs QP');
    ylim([0 100]);  % Set y-axis limits from 0 to 100

    % Plot percentage of splits vs Bitrate
    subplot(1, 2, 2);
    plot(bits_per_frame, split_percentages, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('Bitrate (bits/frame)');
    ylabel('Percentage of Split Blocks (%)');
    title('Block Split Decisions vs Bitrate');
    ylim([0 100]);  % Set y-axis limits from 0 to 100

    % Save the plots
    saveas(gcf, '../Outputs/vbs_statistics.png');

    % Print statistics
    fprintf('\nVBS Statistics:\n');
    fprintf('==============\n');
    for i = 1:length(QPs)
        fprintf('QP = %d:\n', QPs(i));
        fprintf('  Split blocks: %.2f%%\n', split_percentages(i));
        fprintf('  Bits per frame: %.2f\n', bits_per_frame(i));
        fprintf('\n');
    end

    % Save numerical results
    save('../Outputs/vbs_statistics.mat', 'QPs', 'split_percentages', 'bits_per_frame');
end

% Add the exp_golomb_decode function here
function data = exp_golomb_decode(encoded)
    % Decode a 1D array of binary values (0s and 1s) encoded using Exponential-Golomb coding
    data = [];  % Initialize an empty array to store decoded values
    idx = 1;
    
    while idx <= length(encoded)
        % Count leading zeros
        leading_zeros = 0;
        while idx <= length(encoded) && encoded(idx) == 0
            leading_zeros = leading_zeros + 1;
            idx = idx + 1;
        end
        
        % Read remainder of binary value
        if idx <= length(encoded) && encoded(idx) == 1
            idx = idx + 1;
            bin_code_length = leading_zeros;
            bin_code = [1, encoded(idx : idx + bin_code_length - 1)];
            idx = idx + bin_code_length;
            
            % Convert to decimal
            value = bin2dec(char(bin_code + '0'));
            
            % Map back to original value
            if mod(value, 2) == 0
                x = (value / 2);
            else
                x = -(value - 1) / 2;
            end
            
            data = [data, x];
        end
    end
end