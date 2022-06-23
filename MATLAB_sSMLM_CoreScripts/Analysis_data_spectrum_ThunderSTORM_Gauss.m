%% Wrapper for spectral_pairFinding based on ThunderSTORM-gaussian input
%Default values:
% angrange = [-0.04*pi 0*pi];
% distrange = [2500 4500];

function [particles, allpossibilities] = Analysis_data_spectrum_ThunderSTORM_Gauss(inputThSTORMcsv,angrange,distrange,savecsv,csvloc)
ThSTORMdata = csvread(inputThSTORMcsv,1,0);
frxylist = ThSTORMdata(:,[2 3 4 5 6 7 8 9 11]); %Gauss THSTORM
%Frxylist is now frame	x [nm]	y [nm]	z [nm]	sigma1 [nm]	sigma2 [nm]	intensity [photon]	offset [photon]	bkgstd [photon]	uncertainty_xy [nm]

[particles, allpossibilities] = spectral_pairFinding_Gauss(frxylist,angrange,distrange);


%%
if savecsv
    ThSTORMdatafull = [particles];

    array = zeros(size(ThSTORMdatafull,1),6);

    array(:,1) = [1:size(ThSTORMdatafull,1)]'; % id
    %Check for 3d
    if median(ThSTORMdatafull(:,[4])) ~= 0
        array(:,[2 3 4 5 6 7]) = ThSTORMdatafull(:,[1 2 3 4 5 6]); %frame x y ang dist 
        %export to OUT file in same folder, same name
        temp = num2cell('"id","frame","x [nm]","y [nm]","z [nm]","ang","dist"');
        dlmwrite(csvloc,temp, '')
        dlmwrite(csvloc,array,'precision',10,'-append');
    else
        array(:,[2 3 4 5 6]) = ThSTORMdatafull(:,[1 2 3 5 6]); %frame x y ang dist 

        %export to OUT file in same folder, same name
        temp = num2cell('"id","frame","x [nm]","y [nm]","ang","dist"');
        dlmwrite(csvloc,temp, '')
        dlmwrite(csvloc,array,'precision',10,'-append');
    end
end
end
