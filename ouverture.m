function IM = ouverture(img)
    imE = erosion(img);
    IM = dilatation(imE);
end