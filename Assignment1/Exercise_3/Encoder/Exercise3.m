addpath('../Utils');  % For utils functions
% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                 % Number of frames to process
blockSize = 8;                   % Block size for motion estimation
searchRange = 8;                 % Search range r = 1,4, and 8



% Pre-process

dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoderEx3(paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange)

decoderEx3(numFrames, paddedWidth, paddedHeight, blockSize, searchRange)
%decoder
compareYUVFrames('../Outputs/Y_only_foreman.yuv', '../Outputs/decoded_Y_foreman.yuv', 352, 288, 10);