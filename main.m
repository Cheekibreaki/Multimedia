addpath('../Utils');  % For utils functions
addpath('../Outputs');  % For utils functions
addpath('../Decoder');  % For Decoder functions
% Parameters
filename = '../foreman_cif-1.yuv';  % YUV file to read
outputFile = '../Outputs/Y_only_foreman.yuv'; % File to store Y-only components
paddedOutputFile = '../Outputs/padded_Y_foreman.yuv'; % File to store padded Y components
referenceFile = '../Outputs/referenceFrames.yuv';
decodedFile = '../Outputs/decoded_Y_foreman.yuv';
width = 352;                     % Frame width
height = 288;                    % Frame height
numFrames = 3;                 % Number of frames to process
searchRange = 8;                 % Search range r = 1,4, and 8
QP = 3;
j = 3;
VBSEnable = true;
FMEEnable = false;
FastME = true;

if(VBSEnable == true)
    j = j-1;
end

blockSize = 2^j;                   % Block size for motion estimation
dct_blockSize = 2^j;
I_Period = 10; 
nRefFrames = 4;                 % Can take value from 1 to 4

lambda = 1;
% Pre-process



dumpYComponentsToFile(filename, width, height, numFrames, outputFile);

[paddedWidth,paddedHeight] = padYComponentsFromFile(outputFile, numFrames, width, height, blockSize, paddedOutputFile);

% encoder
encoder(referenceFile, paddedOutputFile, numFrames,paddedWidth, paddedHeight, blockSize, searchRange, dct_blockSize, QP, I_Period, nRefFrames,lambda,VBSEnable, FMEEnable,FastME );
[total_byte,bytes_list] = decoder(decodedFile);
%decoder
compareYUVFrames(referenceFile, outputFile, decodedFile, width, height, numFrames);
calculatePSNR(decodedFile, paddedOutputFile, width, height, numFrames)



%All the graph for report:
%generate_rd_analysis();

analyze_vbs_statistics();

%analyze_reference_frames();

%test_lambda_qp()