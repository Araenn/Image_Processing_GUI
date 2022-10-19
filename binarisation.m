function IM = binarisation(img, val)
    img = im2gray(img);
    [nl, nc] = size(img);
    for i = 1:nl
        for j = 1:nc
            if img(i, j) > val
                IM(i, j) = 1;
            else
                IM(i, j) = 0;
            end
        end
    end
    IM = IM.*255;
end