function IM = erosion(img)
    img = im2gray(img);
    [nc, nl] = size(img);
    lignesSuppl = uint8(255.*ones(1, nl));
    IML = [lignesSuppl;img;lignesSuppl]; %rajout ligne supp
    [nc, nl] = size(IML);
    colSuppl = uint8(255.*ones(nc, 1)); %rajout col supp
    IMC = [colSuppl IML colSuppl];
    IM = IMC;
    [nc, nl] = size(img);
    
    for i = 1:nc
        for j = 1:nl
            pix = [IMC(i, j + 1), IMC(i+1, j), IMC(i+1, j+1), IMC(i+1,j+2), IMC(i+2,j+1)];
            mpix = min(pix(:));
            IM(i, j) = mpix;
        end
    end
    IM(1,:,:) = [];
    IM(end,:,:) = [];
    IM(:,1,:) = [];
    IM(:,end,:) = [];
    IM = uint8(IM);
end