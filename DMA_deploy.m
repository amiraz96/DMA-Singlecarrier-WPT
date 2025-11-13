function [DMA_element_loc, H, RFC_num, passive_num] = DMA_deploy(f, antenna_length, start_point)

%% Initial Parameters
c0 = physconst('LightSpeed');
lambda = c0/f;

%% Deploy the DMA and Calculate element Loss !

DMA_WG_dist = lambda/2; DMA_element_dist = lambda/5;

RFC_num = floor(antenna_length / (DMA_WG_dist)); 
passive_num = floor(antenna_length / (DMA_element_dist));

H = zeros(passive_num*RFC_num, passive_num*RFC_num);
for i = 1:RFC_num
    for l = 1:passive_num
        H((i-1)*passive_num + l, (i-1)*passive_num + l) = H_DMA(f, ...
            antenna_length, (l - 1)*DMA_element_dist);
    end
end

DMA_loc = [start_point(1) - floor(passive_num/2)*DMA_element_dist ...
    start_point(2) - floor(RFC_num/2)*DMA_WG_dist start_point(2)];

DMA_element_loc = zeros(RFC_num*passive_num, 3);

for i = 1:RFC_num
    for j = 1:passive_num
        DMA_element_loc((i-1)*passive_num + j, :) = [DMA_loc(1) + (j - 1)*DMA_element_dist ...
            DMA_loc(2) + (i - 1)*DMA_WG_dist, DMA_loc(3)];
    end
end

end