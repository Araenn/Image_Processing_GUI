function IM = img_ouverture()
    name_file=uigetfile( '*.')
    Fid = fopen(name_file);
    L = inputdlg('Lignes');
    C = inputdlg('Colonnes');
    B = inputdlg('Bandes');
    L = str2double(L);
    B = str2double(B);
    C = str2double(C);
    donnees = fread(Fid, C * L * B,'*ubit8');
    fclose(Fid);
    donnees = reshape(donnees, C, L, B);
    donnees(:,:,1)=donnees(:,:,1)';
    donnees(:,:,2)=donnees(:,:,2)';
    donnees(:,:,3)=donnees(:,:,3)';
    IM(:,:,1)=donnees(:,:,3);
    IM(:,:,2)=donnees(:,:,2);
    IM(:,:,3)=donnees(:,:,1);
    imshow(IM)
end