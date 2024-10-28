addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                 % Number of frames to process
blockSize = 8;                   % Block size for motion estimation
searchRange = 2;                 % Search range r = 1,4, and 8
dct_blockSize = 2;
QP = 3;
I_Period = 1; 


% Pre-process

dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoderEx4(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period)

[total_byte,bytes_list] = decoderEx4(decodedFile);
%decoder
%compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);
calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames)

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

