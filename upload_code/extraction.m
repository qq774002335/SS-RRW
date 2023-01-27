function [extracted_watermark,recover_image] = extract(attacked_image,p,g1)
    %EXTRACT 此处显示有关此函数的摘要
    %% EXTRACTION
    % 1. DWT
    % attacked_image = double(attacked_image);
    [LLw HLw LHw HHw] = dwt2(attacked_image, 'haar');

    % 2.Blocks
    [LH_row LH_col] = size(LHw);
    row_blk_num = LH_row/8;
    col_blk_num = LH_col/8;
    blockw = zeros(8,8);
    Aw = zeros(64,64);
    blocks_num = 0;
    for i = 1:row_blk_num
        for j = 1:col_blk_num
            blocks_num = blocks_num+1;
            blockw = LHw((i - 1) * 8 + 1 : i * 8, (j - 1) * 8 + 1 : j * 8);
            for m = 1:8
                Aw(blocks_num,(m-1)*8+1:m*8) = blockw(m,:);
            end

        end
    end
    %Extract random p vectors from the array of positions
    array11_get = zeros(1024,16);
    for n = 1:blocks_num
        for i = 1:16
            array11_get(n,i) = Aw(n,p(n,i));
        end
    end

    % 3. Extract parameter
    N_data = load('N_data.mat');
    N0 = N_data.N0;
    N1 = N_data.N1;
    T_data = load('length.mat');
    Tmax=T_data.Tmax;
    G=T_data.G;
    
    % 4. Recognize watermark value
    array_extract = zeros(1024,16);
    watermark_extract = zeros(1,1);
    Z1 = zeros(1,1024);
    watermark_col=32;
    watermark_row=32;
    count = 0;
    count_range=0;
    for i = 1:watermark_col
        for j = 1:watermark_row
            count = count+1; 
            Z1(count) = array11_get(count,:)*g1'/16;
            if Z1(count)<=(2*Tmax+G)/2 &&  Z1(count)>= -(2*Tmax+G)/2
                count_range = count_range+1;
            end
        end
    end

    T = (g1*g1')/16;


    if count_range == N0  
        %None/Minor attack
        count = 0;
        for i = 1:watermark_row
            for j = 1:watermark_col
                count = count+1; 
                if ((array11_get(count,:)*g1(1,:)')/16 >= (2*Tmax+G)/2)
                    watermark_extract((i-1)*32+j) = 1;
                    array_extract(count,:) = array11_get(count,:) - floor((Tmax+G)/T*g1);
                elseif((array11_get(count,:)*g1(1,:)')/16  <= -(2*Tmax+G)/2)
                    watermark_extract((i-1)*32+j) = 1;
                    array_extract(count,:) = array11_get(count,:) + floor((Tmax+G)/T*g1);
                else
                    watermark_extract((i-1)*32+j) = 0;
                    array_extract(count,:) = array11_get(count,:);
                end
            end
        end

    else
        %Severe attack
        B1(:) = sort(abs(Z1(:)));
        lim_0 = B1(N0);
        lim_1 = B1(N0+1);
        T1 = lim_0;
        G = lim_1 - lim_0;
        count = 0;
        for i = 1:watermark_row
            for j = 1:watermark_col
                count = count+1; 
                if ((array11_get(count,:)*g1(1,:)')/16 >= lim_1)
                    watermark_extract((i-1)*32+j) = 1;
                    array_extract(count,:) = array11_get(count,:) - floor((T1+G)/T*g1);
                elseif((array11_get(count,:)*g1(1,:)')/16  <= -lim_1)
                    watermark_extract((i-1)*32+j) = 1;
                    array_extract(count,:) = array11_get(count,:) + floor((T1+G)/T*g1);
                else
                    watermark_extract((i-1)*32+j) = 0;
                    array_extract(count,:) = array11_get(count,:);
                end
            end
        end
    end

    % 5. Reversibly restore pixel values
    A_extract = Aw;
    for i = 1:1024
        for j = 1:16
            if array_extract(i,j) ~= 0
                A_extract(i,p(i,j)) = array_extract(i,j);
            end
        end
    end
    %Restore Image
    [A_row A_col] = size(A_extract);
    new_LH = zeros(256,256);
    count = 0;
    for m = 1:32
        for i = 1:32
            count=count+1;
            for j = 1:8
                new_LH((m-1)*8+j,(i-1)*8+1:i*8)= A_extract(count,(j-1)*8+1:j*8);
            end

        end
    end

    % 6. Display Image
    recover_image = idwt2(LLw, HLw, new_LH, HHw, 'haar');
    recover_image = uint8(recover_image);
    extracted_watermark =watermark_extract;
end

