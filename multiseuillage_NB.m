function IM = multiseuillage_NB(img, nNB)
    img = im2gray(img);
    [nl, nc] = size(img);
    nPixelSeuil = 255/nNB;
    newVal = zeros(nNB, 1);
    newVal(1) = 0;
    newVal(nNB) = 255;
    for i = 2:nNB
        newVal(i) = 255/(nNB*i-1);
    end

    for i = 1:nl
        for j = 1:nc
            if 0 < img(i, j) && img(i, j) < nPixelSeuil
                IM(i, j) = newVal(1);
            end
            for k = 1:nNB
                if nPixelSeuil*k < img(i, j) && img(i, j) < nPixelSeuil*2*k
                    IM(i, j) = newVal(k+1) + (nPixelSeuil*k - newVal(k+1));
                end
            end
        end
    end
    IM = round(IM);
end