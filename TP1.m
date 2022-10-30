clc; clear; close all

% IM(:,:,3) = [11 4 5 4 5 7 4 5;10 6 8 4 5 5 6 7;12 5 7 10 4 6 7 7;
%     13 50 50 50 50 10 13 12;10 50 52 52 51 10 11 12;11 50 51 52 52 10 12 12;
%    12 50 50 50 51 10 11 13;10 10 10 11 12 10 12 13];

IM = imread("jinx_fire.jpg");
imR = IM(:, :, 1);
imG = IM(:, :, 2);
imB = IM(:, :, 3);

figure(1)
imshow(IM)

figure(2)
h1 = histogramme(imR);
h2 = histogramme(imG);
h3 = histogramme(imB);
stem(h1, '.', 'r')
hold on
stem(h2, '.', 'g')
stem(h3, '.', 'b')
grid()

figure(3)
IM2 = IM;
IM2(:, :, 1) = transf_lineair(imR);
IM2(:, :, 2) = transf_lineair(imG);
IM2(:, :, 3) = transf_lineair(imB);
imshow(IM2)

% % figure(4)
% h1 = histogramme(IM2(:, :, 1));
% h2 = histogramme(IM2(:, :, 2));
% h3 = histogramme(IM2(:, :, 3));
% % stem(h1, '.', 'r')
% % hold on
% % stem(h2, '.', 'g')
% % stem(h3, '.', 'b')
% % grid()
% 
% % figure(5)
% IM3(:, :, 1) = egalisation(IM2(:, :, 1), h1);
% IM3(:, :, 2) = egalisation(IM2(:, :, 2), h2);
% IM3(:, :, 3) = egalisation(IM2(:, :, 3), h3);
% % imshow(IM3)
% 
% %Prewitt
% 
% [ImMix, IMvert, IMhor] = contours(IM3, 1);
% 
% 
% figure(4)
% imshow(ImMix)
% 
% %Sobel
% [ImMix2, IMvert, IMhor] = contours(IM3, 2);
% 
% figure(5)
% imshow(ImMix2)
% 
% sigma = 3;
% IM5(:, :, 1) = gauss(IM2(:, :, 1), sigma, 5);
% IM5(:, :, 2) = gauss(IM2(:, :, 2), sigma, 5);
% IM5(:, :, 3) = gauss(IM2(:, :, 3), sigma, 5);
% 
% figure(6)
% imshow(IM5)
% 
% figure(7)
% h1 = histogramme(IM5(:, :, 1));
% h2 = histogramme(IM5(:, :, 2));
% h3 = histogramme(IM5(:, :, 3));
% 
% IM6(:, :, 1) = transf_lineair(IM5(:, :, 1));
% IM6(:, :, 2) = transf_lineair(IM5(:, :, 2));
% IM6(:, :, 3) = transf_lineair(IM5(:, :, 3));
% imshow(IM6)
% 
% IM7 = binarisation(IM2, 128);
% hbin = histogramme(IM7);
% figure(8)
% imshow(IM7)
% 
% figure(9)
% stem(hbin, '.')
% grid()
% 
% figure(10)
% IM8 = multiseuillage_NB(IM2, 25);
% imshow(IM8)

% figure(11)
% hbinNB = histogramme(IM8);
% stem(hbinNB, '.')
% grid()

figure(11)
% SE = [0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0];
% imE = erosion(IM2, SE);
% imshow(imE)

% IM = im2gray(IM);
% [nc, nl] = size(IM);
% lignesSuppl = uint8(255.*ones(1, nl));
% IML = [lignesSuppl;IM;lignesSuppl];
% [nc, nl] = size(IML);
% colSuppl = uint8(255.*ones(nc, 1));
% IMC = [colSuppl IML colSuppl];
% 
% IMF = IMC;
% [nc, nl] = size(IM);
% 
% for i = 1:nc
%     for j = 1:nl
%         pix = [IMC(i, j + 1), IMC(i+1, j), IMC(i+1, j+1), IMC(i+1,j+2), IMC(i+2,j+1)];
%         mpix = min(pix(:));
%         IMF(i, j) = mpix;
%     end
% end
% IMF(1,:,:) = [];
% IMF(end,:,:) = [];
% IMF(:,1,:) = [];
% IMF(:,end,:) = [];
IMF = erosion(IM);
imshow(IMF)

figure(12)
SE = [0 1 0; 1 1 1; 0 1 0];
IMG = imerode(im2gray(IM), SE);
imshow(IMG)