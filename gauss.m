function IM = gauss(im, sigma, T)
    for i = 1:T
        for j = 1:T
            fil(i, j) = exp(-(i^2+j^2)/(2*sigma^2));
        end
    end
    fil = (1/(2*pi*sigma)).*fil;
    IM = conv2(im, fil);
    IM = uint8(IM);
end