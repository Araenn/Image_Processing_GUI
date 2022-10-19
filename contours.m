function [IMmix, IMvert, IMhor] = contours(im, c)
    im = rgb2gray(im);
    Hs1 = [1 c 1; 0 0 0; -1 -c -1];
    Hs2 = [1 0 -1; c 0 -c; 1 0 -1];
    IMvert = conv2(im, Hs1);
    IMvert = uint8(IMvert);
    IMhor = conv2(im, Hs2);
    IMhor = uint8(IMhor);
    IMmix = abs(IMvert) + abs(IMhor);
    IMmix = uint8(IMmix);
end