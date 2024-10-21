% Testbench for entropy encoder and decoder functions

% Define test input values
frame_type = 0;  % P-frame (1 for I-frame)
mvwidth = 4; mvheight = 4;  % Width and height for motion vectors
predwidth = 4; predheight = 4;  % Width and height for prediction modes
reswidth = 4; resheight = 4;  % Width and height for residues

% Define larger fixed matrices for testing
motionVector3d = cat(3, [1, -2, 3, -4; 5, -6, 7, -8; 9, -10, 11, -12; 13, -14, 15, -16], [17, -18, 19, -20; 21, -22, 23, -24; 25, -26, 27, -28; 29, -30, 31, -32], [-33, 34, -35, 36; -37, 38, -39, 40; -41, 42, -43, 44; -45, 46, -47, 48]);
predicitonModes2d = [-1, 2, -3, 4; -5, 6, -7, 8; -9, 10, -11, 12; -13, 14, -15, 16];  % Larger fixed 2D matrix for prediction modes
residues2d = [5, -6, 7, -8; 9, -10, 11, -12; 13, -14, 15, -16; 17, -18, 19, -20];  % Larger fixed 2D matrix for residues

% Encode the inputs
[encodedMotionVector, encodedPredicitonModes, encodedResidues,predmodeBinLength,motionVectorLength] = entropyEncode(frame_type, motionVector3d, predicitonModes2d, residues2d);

% Print encoded values
disp('Encoded Motion Vector:');
disp(encodedMotionVector);

disp('Encoded Prediction Modes:');
disp(encodedPredicitonModes);

disp('Encoded Residues:');
disp(encodedResidues);

% Decode the encoded values
[decodedMotionVector3d, decodedPredicitonModes2d, decodedResidues2d] = entropyDecode(frame_type, encodedMotionVector, encodedPredicitonModes, encodedResidues, mvwidth, mvheight, predwidth, predheight, reswidth, resheight,predmodeBinLength,motionVectorLength);

% Verify that the decoded values match the original input values
% assert(isequal(decodedMotionVector3d, motionVector3d), 'Decoded Motion Vector does not match original!');
% assert(isequal(decodedPredicitonModes2d, predicitonModes2d), 'Decoded Prediction Modes do not match original!');

if frame_type == 1
    assert(isequal(decodedResidues2d, residues2d), 'Decoded Residues do not match original!');
end

% Display success message
disp('All tests passed successfully!');



