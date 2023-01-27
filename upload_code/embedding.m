function [watermarked_image,p,g1] = embedding(cover_image,watermark_logo,G)
    %EMBEDDING 
    N1 = 0;N0 =0;
    for i = 1:32
        for j = 1:32
            if watermark_logo((i-1)*32+j) == 1
                N1 = N1+1;
            else
                N0 = N0+1;
            end
        end
    end
    save('N_data.mat','N0','N1');
    %% Arnold Transformation
    % arnold_iter = 1;
    % [ watermark_arnold ] = arnold( watermark_logo, arnold_iter );
    watermark_arnold = watermark_logo;

    %% DWT
    % Apply Haar wavelet and decompose cover image into four sub-bands
    cover_image = double(cover_image);
    [LL, HL, LH, HH] = dwt2(cover_image, 'haar');


    %Depart LH to block
    [y_row y_col] = size(LH);
    row_blk_num = y_row/8;
    col_blk_num = y_col/8;
    data = load('data.mat');
    g1 = data.g1;
    p = data.p;
    A = zeros(1024,64);
    blocks_num = 0;
    for i = 1:row_blk_num
        for j = 1:col_blk_num
            blocks_num = blocks_num+1;
            block = LH((i - 1) * 8 + 1 : i * 8, (j - 1) * 8 + 1 : j * 8);
            for m = 1:8
                A(blocks_num,(m-1)*8+1:m*8) = block(m,:);
            end        
            [B(blocks_num,:),I(blocks_num,:)] = sort(abs(A(blocks_num,:)));
        end
    end

    blocks_num = row_blk_num * col_blk_num; 
    p = zeros(blocks_num,16);
    array1 = zeros(blocks_num,16);

    %According to the position vector p, 16 random vectors specified by each block are arranged and stored in array1
    for n = 1:blocks_num

        row_num = size(A(1,:),2);
        p(n,1:16) = I(n,1:16);
        for i = 1:16
            array1(n,i) = A(n,p(n,i));
        end
    end

    %% Embedded watermarking
    [watermark_row watermark_col] = size(watermark_arnold);
    count = 0;
    array11 = zeros(blocks_num,16);

    T = (g1*g1')/16;
    m=1;n=1;Z=zeros(1,1);ZZ=zeros(1,1);
    for i = 1:watermark_row
        for j = 1:watermark_col
            count = count+1; 
            Z(count) = array1(count,:)*g1'/16;
        end
    end
    Tmax = max(Z);
    Tmin = min(Z);


    count = 0;m=1;n=1;
    positive0 = 0;
    negative0 = 0;
    %Store the maximum and minimum range
    save('length.mat','Tmax','Tmin','G');
    for i = 1:watermark_row
        for j = 1:watermark_col
            count = count+1; 

            if(array1(count,:)*g1'/16 >= 0)
                if watermark_arnold((i-1)*32+j) == 1
                    array11(count,:) = array1(count,:) + floor((Tmax+G)/T.*g1); 
                else
                    array11(count,:) = array1(count,:);
                    positive0 = positive0 + 1;
                end            
            elseif((array1(count,:)*g1'/16 < 0))
                if watermark_arnold((i-1)*32+j) == 1
                   array11(count,:) = array1(count,:) - floor((abs(Tmax)+G)/T.*g1); 
                else
                    array11(count,:) = array1(count,:);
                    negative0 = negative0 + 1;
                end            
            end

            ZZ(count) = array11(count,:)*g1'/16;

        end
    end

    %Restores the modified pixel to the specified position
    for i = 1:1024
        for j = 1:16
            if array11(i,j) ~= 0
                A(i,p(i,j)) = array11(i,j);
            end
        end
    end
    %Restore Image
    [A_row A_col] = size(A);
    new_LH = zeros(256,256);
    count = 0;
    for m = 1:32
        for i = 1:32
            count=count+1;
            for j = 1:8
                new_LH((m-1)*8+j,(i-1)*8+1:i*8)= A(count,(j-1)*8+1:j*8);
            end        
        end
    end

    watermarked_image = idwt2(LL, HL, new_LH, HH, 'haar');
    watermarked_image = uint8(watermarked_image);

end

