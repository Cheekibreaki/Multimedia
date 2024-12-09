function per_block_row_budget = get_row_budget_base_on_frame(per_frame_bit_budget,paddedWidth,paddedHeight,block_width)
    full_pixel_size = paddedWidth * paddedHeight;
    block_size = block_width * block_width;

    per_block_budget = per_frame_bit_budget/(full_pixel_size/block_size);

    per_block_row_budget = per_block_budget * (paddedWidth/block_width);
end