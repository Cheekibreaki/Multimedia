addpath('../Utils');  % For utils functions
addpath('../Outputs');
% Parameters
filename = '../../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
width = 352;
height = 288;
numFrames = 10;  
blockSize = 8; % motion estinmation
dct_blockSize = 8;
searchRange = 8;
QP = 6;
I_Period = 1; 

% Pre-process

dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoderEx4(filename, numFrames, width, height, blockSize, dct_blockSize, searchRange, QP, I_Period);
fprintf('Encoder finished.\n');
%decoder
decoderEx4(numFrames, width, height, blockSize, searchRange, QP, I_Period);
fprintf('Decoder finished.\n');

%compareYUVFrames('../Outputs/Y_only_foreman.yuv', '../Outputs/decoded_Y_foreman.yuv', 352, 288, 10);A

