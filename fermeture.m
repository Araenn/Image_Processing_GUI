function IM = fermeture(img)
    imD = dilatation(img);
    IM = erosion(imD);
end