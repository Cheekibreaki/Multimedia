function Q = qp(i, QP)
    % Check if QP is in the valid range
    if QP < 0 || QP > log2(i) + 7
        error('QP must be between 0 and log2(i) + 7 (inclusive)');
    end
    
    % Initialize the quantization matrix Q
    Q = zeros(i, i);
    for row = 1:i
        for col = 1:i
            Q(row, col) = 2 ^ (QP + floor((row + col - 2) / i));
        end
    end
end

% Example Usage
i1 = 2; QP1 = 0;
fprintf('Q for i=2 and QP=0:\n');
disp(qp(i1, QP1));

fprintf('\n');
i2 = 4; QP2 = 2;
fprintf('Q for i=4 and QP=2:\n');
disp(qp(i2, QP2));
