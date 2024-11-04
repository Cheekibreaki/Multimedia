classdef VBS
    properties
        VBS_matrix  % Matrix storing the chosen block type (1 = single MV, 2 = four sub-block MVs)
        VBSEnable   % Parameter for enabling/disabling VBS
    end
    
    methods
        function obj = Create_VBS_matrix(obj, width, height, j, VBSEnable)
            % Calculate the matrix dimensions based on block size
            block_size = 2^(j-1);
            matrix_width = width / block_size;
            matrix_height = height / block_size;
            
            % Initialize VBS_matrix
            obj.VBS_matrix = ones(matrix_height, matrix_width);  % Default to single MV (1)
            
            if VBSEnable
                % Iterate over each block and perform RDO
                for i = 1:matrix_height
                    for k = 1:matrix_width
                        % Ensure the block position matches the larger 2^j x 2^j requirement
                        if mod(i, 2) == 1 && mod(k, 2) == 1
                            % Calculate RD Cost for single block (2^j x 2^j)
                            single_block_cost = obj.compute_rd_cost();

                            % Calculate RD Cost for four sub-blocks
                            four_sub_blocks_cost = obj.compute_4rd_cost();

                            % Compare costs and set VBS_matrix accordingly
                            if four_sub_blocks_cost < single_block_cost
                                % Set the 4 corresponding sub-blocks to "2" (four separate MVs)
                                obj.VBS_matrix(i, k) = 2;
                                obj.VBS_matrix(i + 1, k) = 2;
                                obj.VBS_matrix(i, k + 1) = 2;
                                obj.VBS_matrix(i + 1, k + 1) = 2;
                            end
                        end
                    end
                end
            else
                % If VBSEnable is false, all blocks will default to 1 (single MV)
                obj.VBS_matrix = ones(matrix_height, matrix_width);
            end
        end
        
        function rd_cost = compute_rd_cost(~)
            % Placeholder function to compute the RD Cost for a single block
            % For demonstration purposes, we're using a random value
            % In real cases, this would involve calculating distortion + rate
            rd_cost = rand() * 100;  % Random value representing RD cost
        end
        
        function rd_cost = compute_4rd_cost(~)
            % Placeholder function to compute the RD Cost for four sub-blocks
            % For demonstration purposes, we're using a random value
            % In real cases, this would involve calculating distortion + rate
            rd_cost = rand() * 100;  % Random value representing RD cost
        end
    end
end

