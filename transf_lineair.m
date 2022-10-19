function IM2 = transf_lineair(img)
    ma = max(max(max(img)));
    ma = double(ma);
    mi = min(min(min(img)));
    mi = double(mi);

    alpha = double((255/(ma-mi)));
    beta = double(((-255*mi)/(ma-mi)));

    IM2 = round(alpha.*img + beta);
end