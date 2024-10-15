function testDCTBlock
    % Input block from the image
    block = [ 5, 11,  8, 10;
              9,  8,  4, 12;
              1, 10, 11,  4;
             19,  6, 15,  7 ];

    % Compute DCT of the block
    dctBlock = dct2(double(block));

    % Expected result (You can compute this with dct2 or verify it manually)
    % Here is a placeholder for the expected DCT result (this will need to be filled in)


end

testDCTBlock
