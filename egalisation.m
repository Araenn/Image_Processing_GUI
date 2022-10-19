function IM2 = egalisation(img, h)
    P = h./sum(h);
    S = zeros(length(P), 1);
    S(1) = P(1);
    ma = max(max(max(img)));
    ma = double(ma);
    
    for i = 2:ma + 1
        S(i) = P(i) + S(i - 1);
    end
    for i = 2:ma
        if P(i - 1) == 0
            S(i - 1) = 0;
        end
    end
    CN = round(S.*255);

%     nnh = zeros(length(P), 1);
%     for i = 2:ma + 1
%         if CN(i) ~= 0
%             nnh(CN(i) + 1) = h(i);
%         end
%     end
    [nl, nc] = size(img);
    IM2 = zeros(nl, nc);

    for i = 1:nl
        for j = 1:nc
            IM2(i, j) = CN(img(i, j) + 1);
        end
    end
    IM2 = uint8(IM2);
end