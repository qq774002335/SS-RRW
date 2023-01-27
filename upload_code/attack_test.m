function [attacked_image, extracted_watermark,recover_image] = attack_test(watermarked_image,attack,param,p,g1)
%TEST3 此处显示有关此函数的摘要
%% Attacks
attacked_image = Attacks(watermarked_image,attack,param);

% %逆旋转
% attack_level = param;
% attacked_image1 = imrotate(attacked_image,-attack_level,'bilinear','crop');
% [f c]=find(attacked_image1);%寻找出所有非零元素的位置,f是横坐标向量,c是纵坐标向量。
% attacked_image1=attacked_image1(min(f):max(f),min(c):max(c));%将图像四周的白点去掉
% [m,n]=size(attacked_image1);
% if(attack_level~=90)
%     attacked_image1=attacked_image1(2:n-2,2:n-2);%去掉第一列和最后一列
% end
% attacked_image = imresize(attacked_image1, [512 512], 'bilinear'); 

% attacked_image = imresize(attacked_image, [512 512], 'bilinear'); 
% psnr_attack = psnr(cover_image,attacked_image);
%% Extract
[extracted_watermark,recover_image] = extract(attacked_image,p,g1);
end

