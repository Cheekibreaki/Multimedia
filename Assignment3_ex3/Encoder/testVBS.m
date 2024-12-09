function test_VBS()
    % Test configuration
    width = 32;
    height = 32;
    j = 4;  % Only one j value for testing

    % Loop through VBSEnable as true and false
    for VBSEnable = [true]
        if VBSEnable
            fprintf('VBSEnable = true\n');
        else
            fprintf('VBSEnable = false\n');
        end

        % Create a VBS object and generate the VBS matrix
        vbs = VBS();
        vbs = vbs.Create_VBS_matrix(width, height, j, VBSEnable);

        % Display the resulting VBS matrix
        fprintf('VBS Matrix for j = %d, VBSEnable = %d:\n', j, VBSEnable);
        disp(vbs.VBS_matrix);

        % Check the dimensions of VBS_matrix
        block_size = 2^(j - 1);
        matrix_width = width / block_size;
        matrix_height = height / block_size;

        assert(size(vbs.VBS_matrix, 1) == matrix_height, 'Incorrect matrix height');
        assert(size(vbs.VBS_matrix, 2) == matrix_width, 'Incorrect matrix width');

        % Verify that the elements of the VBS matrix are either 1 or 2
        assert(all(ismember(vbs.VBS_matrix(:), [1, 2])), 'Matrix should only contain 1s or 2s');

        % If VBSEnable is false, ensure all elements are 1
        if ~VBSEnable
            assert(all(vbs.VBS_matrix(:) == 1), 'All elements should be 1 when VBSEnable is false');
        end

        % Display a message that this test configuration passed
        fprintf('Test passed for j = %d, VBSEnable = %d\n\n', j, VBSEnable);
    end
end

test_VBS