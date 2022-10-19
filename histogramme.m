function h = histogramme(img)
    h = zeros(256, 1);
    [nl, nc] = size(img);
    for i = 1:nl
        for j = 1:nc
            val = img(i, j) + 1;
            h(val) = h(val) + 1;
        end
    end
end