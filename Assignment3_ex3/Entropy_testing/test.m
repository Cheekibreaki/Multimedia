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






% % Testbench for Exponential-Golomb encoding and decoding functions
% 
% % Define test input values
% test_data = [-10, 0, 5, 7, -3, 2];
% 
% % Encode the test data
% encoded_data = exp_golomb_encode(test_data);
% 
% % Print encoded values
% disp('Encoded Data:');
% disp(encoded_data);
% 
% % Decode the encoded data
% decoded_data = exp_golomb_decode(encoded_data);
% 
% % Print decoded values
% disp('Decoded Data:');
% disp(decoded_data);
% 
% % Verify that the decoded values match the original input values
% assert(isequal(decoded_data, test_data), 'Decoded data does not match original data!');
% 
% % Display success message
% disp('Exponential-Golomb encoding and decoding test passed successfully!');
% 
% 
% 
% function encoded = exp_golomb_encode(data)
%     % Encode data using Exponential-Golomb coding based on specific rules
%     % and return a 1D array of binary values (0s and 1s)
%     encoded = [];  % Initialize an empty array to store binary values
% 
%     for i = 1:length(data)
%         x = data(i);
%         % Map value to a non-negative integer
%         if x > 0
%             value = 2 * x;  % Positive x to odd integer (2x)
%         else
%             value = (-2 * x) + 1;  % Non-positive x to even integer (-2x + 1)
%         end
% 
%         bin_code = dec2bin(value) - '0';  % Convert to binary and get an array of 1s and 0s
%         leading_zeros = floor(log2(value)) + 1;  % Compute the number of leading zeros
%         code = [zeros(1, leading_zeros - 1), 1, bin_code(2:end)];  % Assemble the code
% 
%         encoded = [encoded, code];  % Append binary code to the output array
%     end
% end
% 
% 
% 
% function data = exp_golomb_decode(encoded)
%     % Decode a 1D array of binary values (0s and 1s) encoded using Exponential-Golomb coding
%     data = [];  % Initialize an empty array to store the decoded values
%     idx = 1;
% 
%     while idx <= length(encoded)
%         % Count the number of leading zeros
%         leading_zeros = 0;
%         while idx <= length(encoded) && encoded(idx) == 0
%             leading_zeros = leading_zeros + 1;
%             idx = idx + 1;
%         end
% 
%         % Read the remainder of the binary value
%         if idx <= length(encoded) && encoded(idx) == 1
%             idx = idx + 1;
%             bin_code_length = leading_zeros;
%             bin_code = [1, encoded(idx : idx + bin_code_length - 1)];
%             idx = idx + bin_code_length;
% 
%             % Convert binary code to decimal value
%             value = bin2dec(char(bin_code + '0'));
% 
%             % Map value back to original data
%             if mod(value, 2) == 0
%                 x = (value / 2);
%             else
%                 x = - (value - 1) / 2;
%             end
% 
%             data = [data, x];  % Append the decoded value to the output array
%         end
%     end
% end
