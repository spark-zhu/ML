        % Gaussian_Fitting_MASKOx=12.5;
       
         Oy=12.5;
        gap=5;
        Gw=25;
        Gh=25;
        numGx = 8;
        numGy = 4;
        Gx=[Ox:(gap+Gw):Ox+(numGx-1)*(gap+Gw)];
        Gy=flip([Oy:(gap+Gh):Oy+(numGy-1)*(gap+Gh)]);
        [X,Y] = meshgrid(Gx,Gy);
        load('C:\Users\xie.lab.ws1\Documents\MATLAB\ML_code\oversampling\Ch_Map_20161207_right_New.mat')
        for row =1:4
            for col =1:8
            Z(row,col) =   P2PForFeature{968}(Ch_Map_new(row,col)-15);
            end
        end
        Intensity
        createFit_2D(X,Y,Z)
       
        colorbar