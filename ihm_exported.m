classdef ihm_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        IHMappUIFigure                matlab.ui.Figure
        WebcamButton                  matlab.ui.control.Button
        MorphologiemathmatiquePanel   matlab.ui.container.Panel
        FermetureButton               matlab.ui.control.Button
        OuvertureButton               matlab.ui.control.Button
        DilatationButton              matlab.ui.control.Button
        ErosionButton                 matlab.ui.control.Button
        SeuillagePanel                matlab.ui.container.Panel
        Slider_mRGB                   matlab.ui.control.Slider
        Slider_mNB                    matlab.ui.control.Slider
        Slider                        matlab.ui.control.Slider
        MultiseuilsRGBButton          matlab.ui.control.Button
        MultiseuilsNBButton           matlab.ui.control.Button
        BinarisationButton            matlab.ui.control.Button
        SaveButton                    matlab.ui.control.Button
        AprsButton                    matlab.ui.control.Button
        AvantButton                   matlab.ui.control.Button
        FiltragePanel                 matlab.ui.container.Panel
        MorphologiemathmatiqueButtonFiltrage  matlab.ui.control.Button
        Filtre                        matlab.ui.control.DropDown
        Bruit                         matlab.ui.control.DropDown
        ContoursPanel                 matlab.ui.container.Panel
        MorphologiemathmatiqueButtonContours  matlab.ui.control.Button
        Sobel                         matlab.ui.control.DropDown
        Prewitt                       matlab.ui.control.DropDown
        ContrastesPanel               matlab.ui.container.Panel
        EgalisationButton             matlab.ui.control.Button
        TransformationLineaireButton  matlab.ui.control.Button
        OuverturedimageButton         matlab.ui.control.Button
        hist_nb                       matlab.ui.control.UIAxes
        hist_rgb                      matlab.ui.control.UIAxes
        figure_apres                  matlab.ui.control.UIAxes
        image_originale               matlab.ui.control.UIAxes
    end

    
    methods (Access = private)
 
       
        
        function h = histogramme(~, img)
            h = zeros(256, 1);
            [nl, nc] = size(img);
            for i = 1:nl
                for j = 1:nc
                    val = img(i, j) + 1;
                    h(val) = h(val) + 1;
                end
            end
        end
        
        function IM2 = transf_lineair(~, img)
            ma = max(max(max(img)));
            ma = double(ma);
            mi = min(min(min(img)));
            mi = double(mi);
        
            alpha = double((255/(ma-mi)));
            alpha = double(alpha);
            beta = double(((-255*mi)/(ma-mi)));
            beta = double(beta);
            IM2 = round(alpha.*img + beta);
        end
        
        function IM2 = egalisation(~, img, h)
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
        
            [nl, nc] = size(img);
            IM2 = zeros(nl, nc);
        
            for i = 1:nl
                for j = 1:nc
                    IM2(i, j) = CN(img(i, j) + 1);
                end
            end
            IM2 = uint8(IM2);
        end
        
        function initModificationHistoric(~)
            
            IM = evalin('base', "IM");
            
            modificationHistoric = [];
            MHCursor = 0;
            MHSize = 0;
            [imageSize,lol] = size(IM);
            
            assignin("base", "imageSize", imageSize);
            assignin("base", "modificationHistoric", modificationHistoric)
            assignin("base", "MHCursor", MHCursor)
            assignin("base", "MHSize", MHSize)        
            
        end
        
        function addModification(~, img)
            
            modificationHistoric = evalin('base', "modificationHistoric");
            MHCursor = evalin('base', "MHCursor");
            MHSize = evalin('base', "MHSize");
            
            if MHCursor == MHSize %Need to append image in array
               modificationHistoric =[modificationHistoric;img];
               MHCursor = MHCursor + 1;
               MHSize = MHSize + 1; 
               assignin("base", "modificationHistoric", modificationHistoric)
            else % write in existing cell and delete old/deprecated next modifications
               MHCursor = MHCursor + 1;
               imageSize = evalin("base", "imageSize");
               images = cat(MHSize,modificationHistoric);
               left = (imageSize - 1) * (MHCursor - 1) + 1;
               right = (imageSize - 1) * (MHCursor) + 1;
               images((left:right),:,:) = img;
               imagesReduced = images((1:right),:,:);
               MHSize = MHCursor;
               assignin("base", "modificationHistoric", imagesReduced)
            end
            assignin("base", "MHCursor", MHCursor)
            assignin("base", "MHSize", MHSize)
           
        end
        
        function img = ctrlYMH(~)
            modificationHistoric = evalin('base', "modificationHistoric");
            MHCursor = evalin('base', "MHCursor");
            MHSize = evalin('base', "MHSize");
            imageSize = evalin("base", "imageSize");
          
            images = cat(MHSize,modificationHistoric);
            if (MHCursor < MHSize)
                MHCursor = MHCursor + 1;
            end
            
            % Don't place that in pre-block
            % It's usefull when no modification has been applied
            left = (imageSize - 1) * (MHCursor - 1) + 1;
            right = (imageSize - 1) * (MHCursor) + 1;
            img = images((left:right),:,:);
            assignin("base", "MHCursor", MHCursor)
        end
        
        function img = ctrlZMH(~)
            modificationHistoric = evalin('base', "modificationHistoric");
            MHCursor = evalin('base', "MHCursor");
            MHSize = evalin('base', "MHSize");
            imageSize = evalin("base", "imageSize");
            
            images = cat(MHSize,modificationHistoric);
            if MHCursor == 1
                img = images((1:(imageSize - 1)),:,:);
            else
                MHCursor = MHCursor - 1;
                left = (imageSize - 1) * (MHCursor - 1) + 1;
                right = (imageSize - 1) * (MHCursor) + 1;
                img = images((left:right),:,:);
                assignin("base", "MHCursor", MHCursor)
            end
            
        end
        
        function IM = binarisation(~, img, val)
            [nl, nc] = size(img);
            IM = img;
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
        
        function IM = multiseuillage_NB(~, img, seuil)
            [nl, nc] = size(img);
            seuil = round(seuil);
            
            listPixel = zeros(nl*nc, 1);
            for l = 1:nl
                for c = 1:nc
                    listPixel(nc * (l - 1) + (c - 1) + 1) = img(l,c);
                end
            end
            
            listPixel = sort(listPixel);
            
            IM = img;
            for l = 1:nl
                for c = 1:nc
                    % Algo de recherche par dichotomie O(log_2(nc * ln))
                    mini = 1;
                    maxi = nc * nl;
                    index = 1;
                    founded = false;
                    while ~founded && mini <= maxi
                        middle = round((maxi + mini) / 2);
                        if listPixel(middle) == img(l, c)
                            founded = true;
                            index = middle;
                        elseif img(l, c) > listPixel(middle)
                            mini = middle + 1;
                        else
                            maxi = middle - 1;
                        end
                    end
                    
                    index = floor((index * seuil) / (nl * nc));
                    grayValue = floor((255 * index) / (seuil - 1));
                    IM(l, c) = grayValue;
                end
            end
            IM = uint8(IM);
        end
        
        function IM = erosion(~, img)
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
        
        function IM = dilatation(~, img)
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
                    mpix = max(pix(:));
                    IM(i, j) = mpix;
                end
            end
            IM(1,:,:) = [];
            IM(end,:,:) = [];
            IM(:,1,:) = [];         
            IM(:,end,:) = [];
            IM = uint8(IM);
        end
        
        function IM = filtre_milieu(~, img, T)
            [nl, nc] = size(img);
            IM = img;
          
            for i = 2:nl-1
                minVect = [img(i - 1, 1);img(i - 1, 2);img(i - 1, 3)];
                maxVect = [img(i - 1, 1);img(i - 1, 2);img(i - 1, 3)];
                
                for j = 1:3
                    currentVectMin = minVect(j);
                    currentVectMax = maxVect(j);
                    for k = 2:3
                        currentPixel = img(i + (k - 2), j);
                        if currentPixel < currentVectMin
                            minVect(j) = currentPixel;
                        elseif currentPixel > currentVectMax
                            maxVect(j) = currentPixel;
                        end
                    end
                end
                    
                minValue = min(minVect(1), min(minVect(2), minVect(3)));
                maxValue = max(maxVect(1), max(maxVect(2), maxVect(3)));
                
                    
                IM(i, 2) = (minValue + maxValue) / 2;
                
                for j = 3:nc-1           
                    minValue = img(i - 1, j + 1);
                    maxValue = img(i - 1, j + 1);
                    for k = 2:3
                        currentPixel = img(i + (k - 2), j + 1);
                        if currentPixel < minValue
                            minValue = currentPixel;
                        elseif currentPixel > maxValue
                            maxValue = currentPixel;
                        end
                    end
                    
                    minVect(1) = minVect(2);
                    minVect(2) = minVect(3);
                    minVect(3) = minValue;
                    
                    maxVect(1) = maxVect(2);
                    maxVect(2) = maxVect(3);
                    maxVect(3) = maxValue;
                    
                    minValue = min(minVect);
                    maxValue = max(maxVect);
                    
                    IM(i, j) = (minValue + maxValue) / 2;
              
                    
                end
            end
       
            
%             for i = 2:nl-1
%                for j = 2:nc-1
%                    v = img(i-1:i+1, j-1:j+1);
%                    v = reshape(v, 1, T^2);
%                    v = sort(v);
%                    IM(i, j) = (v(1)+v(9))/2;
%                end
%             end
        end
        
        function IM = filtre_median(~, img, T)
            [nl, nc] = size(img);
            IM = img;
            for i = 2:nl-1
               for j = 2:nc-1
                   v = [img(i-1,j-1);img(i-1,j);img(i-1,j+1);
                        img(i,j-1);img(i,j);img(i,j+1);
                        img(i+1,j-1);img(i+1,j);img(i+1,j+1);];
                   v = sort(v);
                   IM(i, j) = v(5);
               end
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function: OuverturedimageButton, image_originale
        function OuverturedimageButtonPushed(app, event)
            clc; close all
            name_file=uigetfile( '*.*');
            [filepath, name, ext] = fileparts(name_file);
            if ext == ".img"
                Fid = fopen(name_file);
                L = 1000;
                C = 1000;
                B = 3;
%                 L = str2double(L);
%                 B = str2double(B);
%                 C = str2double(C);
                donnees = fread(Fid, C * L * B,'*ubit8');
                fclose(Fid);
                donnees = reshape(donnees, C, L, B);
                donnees(:,:,1)=donnees(:,:,1)';
                donnees(:,:,2)=donnees(:,:,2)';
                donnees(:,:,3)=donnees(:,:,3)';
                IM(:,:,1)=donnees(:,:,3);
                IM(:,:,2)=donnees(:,:,2);
                IM(:,:,3)=donnees(:,:,1);
            else
                IM = imread(name_file);
            end
            
            imshow(IM, 'Parent', app.image_originale)
            assignin("base", "IM", IM)
            
            imgNB = im2gray(IM);
            hNB = histogramme(app,imgNB);
            app.hist_nb.Visible = "on";
            stem(app.hist_nb, hNB, '.')
            
            hR = histogramme(app, IM(:, :, 1));
            hG = histogramme(app, IM(:, :, 2));
            hB = histogramme(app, IM(:, :, 3));
            
            app.hist_rgb.Visible = "on";
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            initModificationHistoric(app);
            addModification(app, IM);
        end

        % Button pushed function: TransformationLineaireButton
        function TransformationLineaireButtonPushed(app, event)
            img = evalin('base', "IM");
            NewImg(:, :, 1) = transf_lineair(app, img(:, :, 1));
            NewImg(:, :, 2) = transf_lineair(app, img(:, :, 2));
            NewImg(:, :, 3) = transf_lineair(app, img(:, :, 3));
            imshow(NewImg, 'Parent', app.figure_apres)
            
            hR = histogramme(app, NewImg(:, :, 1));
            hG = histogramme(app, NewImg(:, :, 2));
            hB = histogramme(app, NewImg(:, :, 3));
            
            app.hist_nb.Visible = "on";
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            
            assignin("base", "IM", NewImg)
            assignin("base", "hR", hR)
            assignin("base", "hG", hG)
            assignin("base", "hB", hB)
            
            addModification(app, NewImg);
        end

        % Button pushed function: EgalisationButton
        function EgalisationButtonPushed(app, event)
            img = evalin("base", "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            IM2(:, :, 1) = egalisation(app, img(:, :, 1), hR);
            IM2(:, :, 2) = egalisation(app, img(:, :, 2), hG);
            IM2(:, :, 3) = egalisation(app, img(:, :, 3), hB);
            
            imshow(IM2, 'Parent', app.figure_apres)
            
            nhR = histogramme(app, IM2(:, :, 1));
            nhG = histogramme(app, IM2(:, :, 2));
            nhB = histogramme(app, IM2(:, :, 3));
            
            stem(app.hist_nb, nhR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, nhG, '.', 'g')
            stem(app.hist_nb, nhB, '.', 'b')
            hold(app.hist_nb, 'off')
            
            assignin("base", "IM2", IM2)
            
            addModification(app, IM2);
        end

        % Value changed function: Bruit
        function BruitValueChanged(app, event)
            value = app.Bruit.Value;
            im = evalin("base", "IM");
            imshow(im, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            if value == "Bruit"
                IM = im;
            elseif value == "Gaussien"
                IM = imnoise(im, "gaussian");
            elseif value == "Poisson"
                IM = imnoise(im, "poisson");
            elseif value == "Poivre et sel"
                IM = imnoise(im, "salt & pepper");
            end
            IM = uint8(IM);
            imshow(IM, 'Parent', app.figure_apres)
            hR = histogramme(app, IM(:, :, 1));
            hG = histogramme(app, IM(:, :, 2));
            hB = histogramme(app, IM(:, :, 3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            assignin("base", "IM3", IM)
            
            addModification(app, IM);
        end

        % Value changed function: Prewitt
        function PrewittValueChanged(app, event)
            value = app.Prewitt.Value;
            im = evalin("base", "IM");
            imshow(im, 'Parent', app.image_originale)
            im = rgb2gray(im);
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            c = 1;
            Hs2 = [1 c 1; 0 0 0; -1 -c -1];
            Hs1 = [1 0 -1; c 0 -c; 1 0 -1];
            IMhor = conv2(im, Hs2);
            IMvert = conv2(im, Hs1);
            if value == "Prewitt"
                imshow(im, 'Parent', app.figure_apres)
            elseif value == "Vertical"
                IMvert = uint8(IMvert);
                imshow(IMvert, 'Parent', app.figure_apres)
                h = histogramme(app, IMvert);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMvert)
            elseif value == "Horizontal"
                IMhor = uint8(IMhor);
                imshow(IMhor, 'Parent', app.figure_apres)
                h = histogramme(app, IMhor);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMhor)
            elseif value == "Mixte"
                IMmix = abs(IMvert) + abs(IMhor);
                IMmix = uint8(IMmix); 
                imshow(IMmix, 'Parent', app.figure_apres)
                h = histogramme(app, IMmix);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMmix)
            end
            
        end

        % Value changed function: Sobel
        function SobelValueChanged(app, event)
            value = app.Sobel.Value;
            im = evalin("base", "IM");
            imshow(im, 'Parent', app.image_originale)
            im = rgb2gray(im);
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            c = 2;
            Hs2 = [1 c 1; 0 0 0; -1 -c -1];
            Hs1 = [1 0 -1; c 0 -c; 1 0 -1];
            IMvert = conv2(im, Hs1);
            IMhor = conv2(im, Hs2);
            
            if value == "Sobel"
                imshow(im, 'Parent', app.figure_apres)
            elseif value == "Vertical"
                IMvert = uint8(IMvert);
                imshow(IMvert, 'Parent', app.figure_apres)
                h = histogramme(app, IMvert);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMvert)
            elseif value == "Horizontal"
                IMhor = uint8(IMhor);
                imshow(IMhor, 'Parent', app.figure_apres)
                h = histogramme(app, IMhor);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMhor)
            elseif value == "Mixte"
                IMmix = abs(IMvert) + abs(IMhor);
                IMmix = uint8(IMmix); 
                imshow(IMmix, 'Parent', app.figure_apres)
                h = histogramme(app, IMmix);
                stem(app.hist_nb, h, '.')
%                 addModification(app, IMmix)
            end
        end

        % Value changed function: Filtre
        function FiltreValueChanged(app, event)
            value = app.Filtre.Value;
            img = evalin("base", "IM3");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            T = 3;
      
            if value == "Filtre"
                IM = img;
            elseif value == "Passe-Bas"
               fil = (1/9)*ones(T, T); 

               IM(:,:,1) = conv2(img(:,:,1), fil);
               IM(:,:,2) = conv2(img(:,:,2), fil);
               IM(:,:,3) = conv2(img(:,:,3), fil);
               IM = uint8(IM);
            elseif value == "Médian"
                IM(:,:,1) = filtre_median(app, img(:,:,1), T);
                IM(:,:,2) = filtre_median(app, img(:,:,2), T);
                IM(:,:,3) = filtre_median(app, img(:,:,3), T);
                IM = uint8(IM);
                addModification(app, IM);
            elseif value == "Milieu"
                IM(:,:,1) = filtre_milieu(app, img(:,:,1), T);
                IM(:,:,2) = filtre_milieu(app, img(:,:,2), T);
                IM(:,:,3) = filtre_milieu(app, img(:,:,3), T);
                IM = uint8(IM);
                addModification(app, IM);
            end
            imshow(IM, 'Parent', app.figure_apres)
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            assignin("base", "IM2", IM)
        end

        % Button pushed function: AvantButton
        function AvantButtonPushed(app, event)
            IM = ctrlZMH(app);
            
            imshow(IM, "Parent", app.figure_apres)
        end

        % Button pushed function: AprsButton
        function AprsButtonPushed(app, event)
            IM = ctrlYMH(app);
            
            imshow(IM, "Parent", app.figure_apres)
        end

        % Button pushed function: SaveButton
        function SaveButtonPushed(app, event)
            imsave(app.figure_apres)
        end

        % Button pushed function: BinarisationButton
        function BinarisationButtonPushed(app, event)
                app.Slider.Visible = "on";
                app.Slider.Enable = "on";
                app.Slider_mNB.Visible = "off";
                app.Slider_mNB.Enable = "off";
                app.Slider_mRGB.Visible = "off";
                app.Slider_mRGB.Enable = "off";
        end

        % Value changed function: Slider
        function SliderValueChanged(app, event)
            value = app.Slider.Value;
            img = evalin("base", "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            seuil = value;
            
            IMB = binarisation(app, img(:,:,1), seuil);
            
            imshow(IMB, "Parent", app.figure_apres)
            
            hbin = histogramme(app, IMB);
            stem(app.hist_nb, hbin, '.')
        end

        % Button pushed function: MultiseuilsNBButton
        function MultiseuilsNBButtonPushed(app, event)
            app.Slider.Visible = "off";
            app.Slider.Enable = "off";
            app.Slider_mNB.Visible = "on";
            app.Slider_mNB.Enable = "on";
            app.Slider_mRGB.Visible = "off";
            app.Slider_mRGB.Enable = "off";
        end

        % Value changed function: Slider_mNB
        function Slider_mNBValueChanged(app, event)
            value = app.Slider_mNB.Value;
            img = evalin("base", "IM");
            img = im2gray(img);
            imshow(img, 'Parent', app.image_originale)
            h = histogramme(app, img);
            stem(app.hist_rgb, h, '.')
            seuil = value;
            
            IMB = multiseuillage_NB(app, img, seuil);
            imshow(IMB, "Parent", app.figure_apres)
            
            hbin = histogramme(app, IMB);
            stem(app.hist_nb, hbin, '.')
        end

        % Button pushed function: MultiseuilsRGBButton
        function MultiseuilsRGBButtonPushed(app, event)
            app.Slider.Visible = "off";
            app.Slider.Enable = "off";
            app.Slider_mNB.Visible = "off";
            app.Slider_mNB.Enable = "off";
            app.Slider_mRGB.Visible = "on";
            app.Slider_mRGB.Enable = "on";
        end

        % Value changed function: Slider_mRGB
        function Slider_mRGBValueChanged(app, event)
            value = app.Slider_mRGB.Value;
            img = evalin("base", "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            seuil = value;
            
            IMB(:,:,1) = multiseuillage_NB(app, img(:,:,1), seuil);
            IMB(:,:,2) = multiseuillage_NB(app, img(:,:,2), seuil);
            IMB(:,:,3) = multiseuillage_NB(app, img(:,:,3), seuil);
            imshow(IMB, "Parent", app.figure_apres)
            
            hbinR = histogramme(app, IMB(:,:,1));
            hbinG = histogramme(app, IMB(:,:,2));
            hbinB = histogramme(app, IMB(:,:,3));
            stem(app.hist_nb, hbinR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hbinG, '.', 'g')
            stem(app.hist_nb, hbinB, '.', 'b')
            hold(app.hist_nb, 'off')
        end

        % Button pushed function: ErosionButton
        function ErosionButtonPushed(app, event)
            img = evalin("base", "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            IM(:,:,1) = erosion(app, img(:,:,1));
            IM(:,:,2) = erosion(app, img(:,:,2));
            IM(:,:,3) = erosion(app, img(:,:,3));
            imshow(IM, 'Parent', app.figure_apres)
            
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            addModification(app, IM)
        end

        % Button pushed function: DilatationButton
        function DilatationButtonPushed(app, event)
            img = evalin("base", "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = evalin("base", "hR");
            hG = evalin("base", "hG");
            hB = evalin("base", "hB");
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            IM(:,:,1) = dilatation(app, img(:,:,1));
            IM(:,:,2) = dilatation(app, img(:,:,2));
            IM(:,:,3) = dilatation(app, img(:,:,3));
            imshow(IM, 'Parent', app.figure_apres)
            
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            addModification(app, IM)
        end

        % Button pushed function: OuvertureButton
        function OuvertureButtonPushed(app, event)
            img = evalin('base', "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = histogramme(app, img(:,:,1));
            hG = histogramme(app, img(:,:,2));
            hB = histogramme(app, img(:,:,3));
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
           
            imE(:,:,1) = erosion(app, img(:,:,1));
            imE(:,:,2) = erosion(app, img(:,:,2));
            imE(:,:,3) = erosion(app, img(:,:,3));
            IM(:,:,1) = dilatation(app, imE(:,:,1));
            IM(:,:,2) = dilatation(app, imE(:,:,2));
            IM(:,:,3) = dilatation(app, imE(:,:,3));
            
            imshow(IM, 'Parent', app.figure_apres)
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            addModification(app, IM)
        end

        % Button pushed function: FermetureButton
        function FermetureButtonPushed(app, event)
            img = evalin('base', "IM");
            imshow(img, 'Parent', app.image_originale)
            hR = histogramme(app, img(:,:,1));
            hG = histogramme(app, img(:,:,2));
            hB = histogramme(app, img(:,:,3));
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
           
            imD(:,:,1) = dilatation(app, img(:,:,1));
            imD(:,:,2) = dilatation(app, img(:,:,2));
            imD(:,:,3) = dilatation(app, img(:,:,3));
            IM(:,:,1) = erosion(app, imD(:,:,1));
            IM(:,:,2) = erosion(app, imD(:,:,2));
            IM(:,:,3) = erosion(app, imD(:,:,3));
            
            imshow(IM, 'Parent', app.figure_apres)
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            addModification(app, IM)
        end

        % Button pushed function: 
        % MorphologiemathmatiqueButtonFiltrage
        function MorphologiemathmatiqueButtonFiltragePushed(app, event)
            img = evalin('base', "IM3");
            imshow(img, 'Parent', app.image_originale)
            hR = histogramme(app, img(:,:,1));
            hG = histogramme(app, img(:,:,2));
            hB = histogramme(app, img(:,:,3));
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
           
            imE(:,:,1) = erosion(app, img(:,:,1));
            imE(:,:,2) = erosion(app, img(:,:,2));
            imE(:,:,3) = erosion(app, img(:,:,3));
            IM(:,:,1) = dilatation(app, imE(:,:,1));
            IM(:,:,2) = dilatation(app, imE(:,:,2));
            IM(:,:,3) = dilatation(app, imE(:,:,3));
            
            imshow(IM, 'Parent', app.figure_apres)
            hR = histogramme(app, IM(:,:,1));
            hG = histogramme(app, IM(:,:,2));
            hB = histogramme(app, IM(:,:,3));
            stem(app.hist_nb, hR, '.', 'r')
            hold(app.hist_nb, 'on')
            stem(app.hist_nb, hG, '.', 'g')
            stem(app.hist_nb, hB, '.', 'b')
            hold(app.hist_nb, 'off')
            addModification(app, IM)
        end

        % Button pushed function: 
        % MorphologiemathmatiqueButtonContours
        function MorphologiemathmatiqueButtonContoursPushed(app, event)
            img = evalin('base', "IM");
            img = im2gray(img);
            imshow(img, 'Parent', app.image_originale)
            h = histogramme(app, img);
            stem(app.hist_rgb, h, '.')
            
            imE = erosion(app, img);
            imD = dilatation(app, img);
            
            IM = abs(imD - imE);
            imshow(IM, 'Parent', app.figure_apres)
            hE = histogramme(app, IM);
            stem(app.hist_nb, hE, '.')
            addModification(app, IM)
        end

        % Button pushed function: WebcamButton
        function WebcamButtonPushed(app, event)
            vid = videoinput('winvideo','HP truevision hd camera','MJPG_640x480');
%             vidRes = vid.VideoResolution; 
%             nBands = vid.NumberOfBands; 
%             hImage = image(app.image_originale, zeros(vidRes(2), vidRes(1), nBands) ); 
            start(vid)
            IM = getsnapshot(vid);
            stop(vid)
            IM = uint8(IM);
            
            imshow(IM, 'Parent', app.image_originale)
            imgNB = im2gray(IM);
            hNB = histogramme(app,imgNB);
            app.hist_nb.Visible = "on";
            stem(app.hist_nb, hNB, '.')
            
            hR = histogramme(app, IM(:, :, 1));
            hG = histogramme(app, IM(:, :, 2));
            hB = histogramme(app, IM(:, :, 3));
            
            app.hist_rgb.Visible = "on";
            stem(app.hist_rgb, hR, '.', 'r')
            hold(app.hist_rgb, 'on')
            stem(app.hist_rgb, hG, '.', 'g')
            stem(app.hist_rgb, hB, '.', 'b')
            hold(app.hist_rgb, 'off')
            
            initModificationHistoric(app);
            addModification(app, IM);
            assignin("base", "IM", IM)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create IHMappUIFigure and hide until all components are created
            app.IHMappUIFigure = uifigure('Visible', 'off');
            app.IHMappUIFigure.Position = [100 100 953 649];
            app.IHMappUIFigure.Name = 'IHM app';

            % Create image_originale
            app.image_originale = uiaxes(app.IHMappUIFigure);
            app.image_originale.XTick = [];
            app.image_originale.YTick = [];
            app.image_originale.ZTick = [];
            app.image_originale.GridAlpha = 0;
            app.image_originale.Visible = 'off';
            app.image_originale.HandleVisibility = 'callback';
            app.image_originale.ButtonDownFcn = createCallbackFcn(app, @OuverturedimageButtonPushed, true);
            app.image_originale.Position = [24 373 209 195];

            % Create figure_apres
            app.figure_apres = uiaxes(app.IHMappUIFigure);
            app.figure_apres.XTick = [];
            app.figure_apres.YTick = [];
            app.figure_apres.ZTick = [];
            app.figure_apres.GridAlpha = 0.15;
            app.figure_apres.Visible = 'off';
            app.figure_apres.HandleVisibility = 'callback';
            app.figure_apres.Position = [259 373 211 195];

            % Create hist_rgb
            app.hist_rgb = uiaxes(app.IHMappUIFigure);
            app.hist_rgb.XLim = [0 256];
            app.hist_rgb.XColor = [0 0 0];
            app.hist_rgb.XTick = [0 255];
            app.hist_rgb.YColor = [0 0 0];
            app.hist_rgb.ZTick = [];
            app.hist_rgb.Color = 'none';
            app.hist_rgb.TitleFontWeight = 'normal';
            app.hist_rgb.Visible = 'off';
            app.hist_rgb.Position = [34 195 191 148];

            % Create hist_nb
            app.hist_nb = uiaxes(app.IHMappUIFigure);
            app.hist_nb.XLim = [0 256];
            app.hist_nb.XColor = [0 0 0];
            app.hist_nb.XTick = [0 255];
            app.hist_nb.YColor = [0 0 0];
            app.hist_nb.Color = 'none';
            app.hist_nb.TitleFontWeight = 'normal';
            app.hist_nb.GridColor = 'none';
            app.hist_nb.Visible = 'off';
            app.hist_nb.Position = [271 195 191 148];

            % Create OuverturedimageButton
            app.OuverturedimageButton = uibutton(app.IHMappUIFigure, 'push');
            app.OuverturedimageButton.ButtonPushedFcn = createCallbackFcn(app, @OuverturedimageButtonPushed, true);
            app.OuverturedimageButton.Position = [66 583 127 38];
            app.OuverturedimageButton.Text = 'Ouverture d''image';

            % Create ContrastesPanel
            app.ContrastesPanel = uipanel(app.IHMappUIFigure);
            app.ContrastesPanel.Title = 'Contrastes';
            app.ContrastesPanel.Position = [572 479 168 117];

            % Create TransformationLineaireButton
            app.TransformationLineaireButton = uibutton(app.ContrastesPanel, 'push');
            app.TransformationLineaireButton.ButtonPushedFcn = createCallbackFcn(app, @TransformationLineaireButtonPushed, true);
            app.TransformationLineaireButton.Position = [4 59 160 31];
            app.TransformationLineaireButton.Text = 'Transformation Lineaire';

            % Create EgalisationButton
            app.EgalisationButton = uibutton(app.ContrastesPanel, 'push');
            app.EgalisationButton.ButtonPushedFcn = createCallbackFcn(app, @EgalisationButtonPushed, true);
            app.EgalisationButton.Position = [4 23 160 26];
            app.EgalisationButton.Text = 'Egalisation';

            % Create ContoursPanel
            app.ContoursPanel = uipanel(app.IHMappUIFigure);
            app.ContoursPanel.Title = 'Contours';
            app.ContoursPanel.Position = [571 325 169 139];

            % Create Prewitt
            app.Prewitt = uidropdown(app.ContoursPanel);
            app.Prewitt.Items = {'Prewitt', 'Horizontal', 'Vertical', 'Mixte'};
            app.Prewitt.ValueChangedFcn = createCallbackFcn(app, @PrewittValueChanged, true);
            app.Prewitt.Position = [8 83 144 27];
            app.Prewitt.Value = 'Prewitt';

            % Create Sobel
            app.Sobel = uidropdown(app.ContoursPanel);
            app.Sobel.Items = {'Sobel', 'Horizontal', 'Vertical', 'Mixte'};
            app.Sobel.ValueChangedFcn = createCallbackFcn(app, @SobelValueChanged, true);
            app.Sobel.Position = [8 46 144 27];
            app.Sobel.Value = 'Sobel';

            % Create MorphologiemathmatiqueButtonContours
            app.MorphologiemathmatiqueButtonContours = uibutton(app.ContoursPanel, 'push');
            app.MorphologiemathmatiqueButtonContours.ButtonPushedFcn = createCallbackFcn(app, @MorphologiemathmatiqueButtonContoursPushed, true);
            app.MorphologiemathmatiqueButtonContours.Position = [4 8 161 28];
            app.MorphologiemathmatiqueButtonContours.Text = 'Morphologie mathématique';

            % Create FiltragePanel
            app.FiltragePanel = uipanel(app.IHMappUIFigure);
            app.FiltragePanel.Title = 'Filtrage';
            app.FiltragePanel.Position = [571 167 169 140];

            % Create Bruit
            app.Bruit = uidropdown(app.FiltragePanel);
            app.Bruit.Items = {'Bruit', 'Gaussien', 'Poisson', 'Poivre et sel'};
            app.Bruit.ValueChangedFcn = createCallbackFcn(app, @BruitValueChanged, true);
            app.Bruit.Position = [12 81 144 27];
            app.Bruit.Value = 'Bruit';

            % Create Filtre
            app.Filtre = uidropdown(app.FiltragePanel);
            app.Filtre.Items = {'Filtre', 'Passe-Bas', 'Médian', 'Milieu'};
            app.Filtre.ValueChangedFcn = createCallbackFcn(app, @FiltreValueChanged, true);
            app.Filtre.Position = [12 44 144 27];
            app.Filtre.Value = 'Filtre';

            % Create MorphologiemathmatiqueButtonFiltrage
            app.MorphologiemathmatiqueButtonFiltrage = uibutton(app.FiltragePanel, 'push');
            app.MorphologiemathmatiqueButtonFiltrage.ButtonPushedFcn = createCallbackFcn(app, @MorphologiemathmatiqueButtonFiltragePushed, true);
            app.MorphologiemathmatiqueButtonFiltrage.Position = [3 9 161 28];
            app.MorphologiemathmatiqueButtonFiltrage.Text = 'Morphologie mathématique';

            % Create AvantButton
            app.AvantButton = uibutton(app.IHMappUIFigure, 'push');
            app.AvantButton.ButtonPushedFcn = createCallbackFcn(app, @AvantButtonPushed, true);
            app.AvantButton.Position = [25 37 69 32];
            app.AvantButton.Text = {'Avant'; ''};

            % Create AprsButton
            app.AprsButton = uibutton(app.IHMappUIFigure, 'push');
            app.AprsButton.ButtonPushedFcn = createCallbackFcn(app, @AprsButtonPushed, true);
            app.AprsButton.Position = [166 37 68 32];
            app.AprsButton.Text = 'Après';

            % Create SaveButton
            app.SaveButton = uibutton(app.IHMappUIFigure, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @SaveButtonPushed, true);
            app.SaveButton.Position = [328 35 80 37];
            app.SaveButton.Text = 'Save';

            % Create SeuillagePanel
            app.SeuillagePanel = uipanel(app.IHMappUIFigure);
            app.SeuillagePanel.Title = 'Seuillage';
            app.SeuillagePanel.Position = [571 10 170 147];

            % Create BinarisationButton
            app.BinarisationButton = uibutton(app.SeuillagePanel, 'push');
            app.BinarisationButton.ButtonPushedFcn = createCallbackFcn(app, @BinarisationButtonPushed, true);
            app.BinarisationButton.Position = [11 99 142 22];
            app.BinarisationButton.Text = 'Binarisation';

            % Create MultiseuilsNBButton
            app.MultiseuilsNBButton = uibutton(app.SeuillagePanel, 'push');
            app.MultiseuilsNBButton.ButtonPushedFcn = createCallbackFcn(app, @MultiseuilsNBButtonPushed, true);
            app.MultiseuilsNBButton.Position = [11 73 145 24];
            app.MultiseuilsNBButton.Text = 'Multiseuils NB';

            % Create MultiseuilsRGBButton
            app.MultiseuilsRGBButton = uibutton(app.SeuillagePanel, 'push');
            app.MultiseuilsRGBButton.ButtonPushedFcn = createCallbackFcn(app, @MultiseuilsRGBButtonPushed, true);
            app.MultiseuilsRGBButton.Position = [12 48 145 24];
            app.MultiseuilsRGBButton.Text = 'Multiseuils RGB';

            % Create Slider
            app.Slider = uislider(app.SeuillagePanel);
            app.Slider.Limits = [0 255];
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.Enable = 'off';
            app.Slider.Visible = 'off';
            app.Slider.Position = [14 34 141 7];

            % Create Slider_mNB
            app.Slider_mNB = uislider(app.SeuillagePanel);
            app.Slider_mNB.Limits = [0 255];
            app.Slider_mNB.ValueChangedFcn = createCallbackFcn(app, @Slider_mNBValueChanged, true);
            app.Slider_mNB.Enable = 'off';
            app.Slider_mNB.Visible = 'off';
            app.Slider_mNB.Position = [14 34 141 7];

            % Create Slider_mRGB
            app.Slider_mRGB = uislider(app.SeuillagePanel);
            app.Slider_mRGB.Limits = [0 255];
            app.Slider_mRGB.ValueChangedFcn = createCallbackFcn(app, @Slider_mRGBValueChanged, true);
            app.Slider_mRGB.Enable = 'off';
            app.Slider_mRGB.Visible = 'off';
            app.Slider_mRGB.Position = [14 34 141 7];

            % Create MorphologiemathmatiquePanel
            app.MorphologiemathmatiquePanel = uipanel(app.IHMappUIFigure);
            app.MorphologiemathmatiquePanel.Title = 'Morphologie mathématique';
            app.MorphologiemathmatiquePanel.Position = [771 205 162 179];

            % Create ErosionButton
            app.ErosionButton = uibutton(app.MorphologiemathmatiquePanel, 'push');
            app.ErosionButton.ButtonPushedFcn = createCallbackFcn(app, @ErosionButtonPushed, true);
            app.ErosionButton.Position = [9 128 142 26];
            app.ErosionButton.Text = 'Erosion';

            % Create DilatationButton
            app.DilatationButton = uibutton(app.MorphologiemathmatiquePanel, 'push');
            app.DilatationButton.ButtonPushedFcn = createCallbackFcn(app, @DilatationButtonPushed, true);
            app.DilatationButton.Position = [9 87 142 28];
            app.DilatationButton.Text = 'Dilatation';

            % Create OuvertureButton
            app.OuvertureButton = uibutton(app.MorphologiemathmatiquePanel, 'push');
            app.OuvertureButton.ButtonPushedFcn = createCallbackFcn(app, @OuvertureButtonPushed, true);
            app.OuvertureButton.Position = [9 51 142 27];
            app.OuvertureButton.Text = 'Ouverture';

            % Create FermetureButton
            app.FermetureButton = uibutton(app.MorphologiemathmatiquePanel, 'push');
            app.FermetureButton.ButtonPushedFcn = createCallbackFcn(app, @FermetureButtonPushed, true);
            app.FermetureButton.Position = [9 11 142 27];
            app.FermetureButton.Text = 'Fermeture';

            % Create WebcamButton
            app.WebcamButton = uibutton(app.IHMappUIFigure, 'push');
            app.WebcamButton.ButtonPushedFcn = createCallbackFcn(app, @WebcamButtonPushed, true);
            app.WebcamButton.Position = [291 587 148 34];
            app.WebcamButton.Text = 'Webcam';

            % Show the figure after all components are created
            app.IHMappUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ihm_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.IHMappUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.IHMappUIFigure)
        end
    end
end