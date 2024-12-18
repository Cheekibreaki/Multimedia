addpath('../Utils');  % For utils functions
addpath('../Decoder');
% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv'
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 10;                  % Number of frames to process
blockSize = 8;                   % Block size for motion estimation
searchRange = 4;                 % Search range r = 1,4, and 8
n = 3;                             % n for rounding factor. 

% Pre-process

dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoderEx3(paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, n)

decoderEx3(decodedFile, numFrames, paddedWidth, paddedHeight, blockSize, searchRange)
%decoder
compareYUVFrames(referenceFile, '../Outputs/Y_only_foreman.yuv', '../Outputs/decoded_Y_foreman.yuv', 352, 288, 10);
