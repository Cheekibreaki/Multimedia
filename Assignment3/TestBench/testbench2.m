addpath('../Utils');  % For utility functions
addpath('../Outputs');  % For output functions
addpath('../Decoder');  % For decoder functions

% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                  % Number of frames to process              % Block size for motion estimation
searchRange = 2;                 % Search range (r = 1, 4, and 8)

blockSize_List = [8, 16];
I_Period_List = [1, 4, 10];
QP_8_List = [3];
QP_16_List = [4];

% Store PSNR results
results = struct();

% Iterate through different configurations
for I_Period = I_Period_List
    for blockSize = blockSize_List
        dct_blockSize = blockSize;
        
        if blockSize == 8
            qpList = QP_8_List;
        elseif blockSize == 16
            qpList = QP_16_List;
        end
        
        total_bits = zeros(1, length(qpList));
        psnrValues = zeros(1, length(qpList));
        encoder_elapsedTimes = zeros(1, length(qpList));
        decoder_elapsedTimes = zeros(1, length(qpList));
        for i = 1:length(qpList)

            dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

            [paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);


            QP = qpList(i);
            % Print selected parameters
            fprintf('Running encoder with parameters: I_Period = %d, blockSize = %d, QP = %d\n', I_Period, blockSize, QP);
            tic;

            % Encode
            encoderEx4(referenceFile, paddedOutputFile, numFrames, width, height, blockSize, searchRange, dct_blockSize, QP, I_Period);
            encoder_elapsedTime = toc;
            tic;

            % Decode
            [total_byte,bytes_list] = decoderEx4(decodedFile);
            bytes_list = bytes_list * 8
            decoder_elapsedTime = toc;
            % Calculate PSNR and store the value
            psnrValues(i) = calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames);
            total_bits(i) = total_byte * 8;
            encoder_elapsedTimes(i) = encoder_elapsedTime;
            decoder_elapsedTimes(i) = decoder_elapsedTime;
        end
        
        % Store results for plotting
        key = sprintf('I_Period_%d_BlockSize_%d', I_Period, blockSize);
        results.(key) = struct('Bytes_list',bytes_list,'QP_List', qpList,'DecoderElapsedTimes', decoder_elapsedTimes, 'EncoderElapsedTimes', encoder_elapsedTimes, 'PSNR_Values', psnrValues, 'total_bits',total_bits);
    end
end

% Plotting PSNR values for blockSize = 8
figure;
hold on;
keys = fieldnames(results);
for i = 1:length(keys)
    key = keys{i};
    if contains(key, 'BlockSize_8')
        qpList = results.(key).QP_List;
        psnrValues = results.(key).PSNR_Values;
        total_bits = results.(key).total_bits;
        bytes_list = results.(key).Bytes_list;
        % Extract I_Period from the key string
        I_Period_str = extractBetween(key, 'I_Period_', '_BlockSize');
        
        % Plot with I_Period as the DisplayName
        plot([1,2,3,4,5,6,7,8,9,10], bytes_list, '-o', 'DisplayName', sprintf('I_Period %s', I_Period_str{1}));
    end
end
xlabel('Frames');
ylabel('Bitcount/Frame');
title('Frames vs Bitcount/Frame for BlockSize = 8');
legend;
grid on;
hold off;

% Plotting PSNR values for blockSize = 16
figure;
hold on;
keys = fieldnames(results);
for i = 1:length(keys)
    key = keys{i};
    if contains(key, 'BlockSize_16')
        qpList = results.(key).QP_List;
        psnrValues = results.(key).PSNR_Values;
        total_bits = results.(key).total_bits;
        bytes_list = results.(key).Bytes_list;
        % Extract I_Period from the key string
        I_Period_str = extractBetween(key, 'I_Period_', '_BlockSize');
        
        % Plot with I_Period as the DisplayName
        plot([1,2,3,4,5,6,7,8,9,10], bytes_list, '-o', 'DisplayName', sprintf('I_Period %s', I_Period_str{1}));
    end
end
xlabel('Frames');
ylabel('Bitcount/Frame');
title('Frames vs Bitcount/Frame for BlockSize = 16');
legend;
grid on;
hold off;




















function avgPsnr = calculatePSNR(decodedFile, originalFile, width, height, numFrames)
    % Parameters:
    %   originalFile - Path to the original Y only file
    %   decodedFile  - Path to the decoded Y only file
    %   width        - Width of the frames
    %   height       - Height of the frames
    %   numFrames    - Number of frames to compare

    % Open the original and decoded files for reading
    fidOriginal = fopen(originalFile, 'r');
    fidDecoded = fopen(decodedFile, 'r');
    
    totalPsnr = 0;
    
    % Iterate through each frame to compare
    for frameIdx = 1:numFrames
        % Read frames
        originalFrame = fread(fidOriginal, [width, height], 'uint8');
        decodedFrame = fread(fidDecoded, [width, height], 'uint8');
        

        % Calculate MSE
        mse = mean((double(originalFrame(:)) - double(decodedFrame(:))).^2);

        % Calculate PSNR
        if mse > 0
            psnrValue = 10 * log10(255^2 / mse);
        else
            psnrValue = Inf; % Perfect match
        end
        % psnrValue = psnr(originalFrame,decodedFrame)
        
        % Accumulate PSNR
        totalPsnr = totalPsnr + psnrValue;
    end
    
    % Average PSNR over all frames
    avgPsnr = totalPsnr / numFrames;
    fprintf('Average PSNR: %.2f dB\n', avgPsnr);
    
    % Close the files
    fclose(fidOriginal);
    fclose(fidDecoded);
end
