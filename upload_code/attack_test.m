function [attacked_image, extracted_watermark,recover_image] = attack_test(watermarked_image,attack,param,p,g1)
%TEST3 �˴���ʾ�йش˺�����ժҪ
%% Attacks
attacked_image = Attacks(watermarked_image,attack,param);

% %����ת
% attack_level = param;
% attacked_image1 = imrotate(attacked_image,-attack_level,'bilinear','crop');
% [f c]=find(attacked_image1);%Ѱ�ҳ����з���Ԫ�ص�λ��,f�Ǻ���������,c��������������
% attacked_image1=attacked_image1(min(f):max(f),min(c):max(c));%��ͼ�����ܵİ׵�ȥ��
% [m,n]=size(attacked_image1);
% if(attack_level~=90)
%     attacked_image1=attacked_image1(2:n-2,2:n-2);%ȥ����һ�к����һ��
% end
% attacked_image = imresize(attacked_image1, [512 512], 'bilinear'); 

% attacked_image = imresize(attacked_image, [512 512], 'bilinear'); 
% psnr_attack = psnr(cover_image,attacked_image);
%% Extract
[extracted_watermark,recover_image] = extract(attacked_image,p,g1);
end

