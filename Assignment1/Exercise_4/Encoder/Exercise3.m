addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv'
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 3;                 % Number of frames to process
blockSize = 8;                   % Block size for motion estimation
searchRange = 8;                 % Search range r = 1,4, and 8
dct_blockSize = 4;
QP = 6;
I_Period = 1; 


% Pre-process

dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoderEx3(paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period)

decoderEx3(decodedFile,numFrames, paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period)
%decoder
compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);