%% Low-dispersion spectral SMLM data analysis
% Basically identical to spectral_pairFinding, but with extra output for
% ThunderSTORM-pSMLM fitting
function [particles, allpossibilities] = spectral_pairFinding_pSMLM(loc_list,curveAngleinput,linking_range)
%%
% Assumed input: loc_list, consisting of frame, x, y columns
%   linking_range: [min_range max_range]
%   curveAngle: [Min and max angle] to look for localization pairs.
% Output: particles array
%     columns: frame, Xpos, Ypos, Zpos, rotation, distance,
%     spatial/spectral PSF size, spatial loc_unc, spatial intensity
%%
%set min and max range for lookup of combinations - this should be based on
%the calibration curve, and maybe even dependant on the angle with the
%other spot as well.
%If it is a 1-by-2 matrix, it just has min and max range.
%If it's anything else, it's a calibration curve.
global verbose
if verbose
disp([char(datetime),'   Starting spectral linking'])
end

    min_range = linking_range(1);
    max_range = linking_range(2);
    min_angle = curveAngleinput(1);
    max_angle = curveAngleinput(2);
    
guess_possible_linking = 0; %0 if only deterministic, 1 if we want to guess unresolvable linkages
abs_error_from_estimated_distance = 1; %was 1; %maximum difference (pixels) between found distance of 2 spots and estimated distance of 2 spots based on their rotation. Unused if just a 1-by-2 matrix is used for min/max linking distance.

particles = zeros(floor(size(loc_list,1)/2),12); %needs to be zeros matrix
curparticle_id = 1;
%% Start analysis
%loop over all frames
%Quick frame index list lookup
frame_index_list = accumarray([loc_list(:,1)], (1:size(loc_list,1))', [], @(idxlist) {idxlist});

for frame = 1:max(loc_list(:,1))
    if mod(frame,10) == 0
        fprintf('Frame %.0f/%.0f\n',frame,max(loc_list(:,1)))
    end
    %Get localizations on this frame
    t = frame_index_list(frame); tt = t{1};
    loc_array = loc_list(tt,2:3); %only X,Y info in matrix
    all_array = loc_list(tt,:); %all info for later finding    
    
    %Add id
    loc_array = [[1:size(loc_array,1)]' loc_array];
    Angles = zeros(size(loc_array,1),size(loc_array,1));
    %Calculate all euclidian distances between all points
    if size(loc_array,1)>35 %If its incredibly large, pdist is faster
        EuclidianDistances = squareform(pdist(loc_array(:,2:3)));
        for ii = 1:size(loc_array,1)
            for j = 1:size(loc_array,1)
                Angles(ii,j) = atan2((loc_array(j,3)-loc_array(ii,3)),((loc_array(j,2)-loc_array(ii,2))));
            end
        end
    else %else 2 forloops is faster for some reason...
        EuclidianDistances = zeros(size(loc_array,1),size(loc_array,1));
        for ii = 1:size(loc_array,1)
            for j = 1:size(loc_array,1)
                EuclidianDistances(ii,j) = sqrt((loc_array(ii,2)-loc_array(j,2))^2+(loc_array(ii,3)-loc_array(j,3))^2);
                Angles(ii,j) = atan2((loc_array(j,3)-loc_array(ii,3)),((loc_array(j,2)-loc_array(ii,2))));
            end
        end
    end
    %Calculate which distances are possible to have a matching pair, currently
    %depending on min_range and max_range and min_angle and max_angle
    possibilities = cell(size(loc_array,1),4);
    Possibilities_size = [];
    if size(loc_array,1) > 1
        for ii = 1:size(loc_array,1)
            possibilities{ii,1} = ii;
            for j = 1:size(loc_array,1)
                if EuclidianDistances(ii,j) > min_range
                    if EuclidianDistances(ii,j) < max_range
                        if Angles(ii,j) > min_angle
                            if Angles(ii,j) < max_angle
                                possibilities{ii,2} = [possibilities{ii,2} j];
                            end
                        end
                    end
                end
            end
            possibilities{ii,3} = size(possibilities{ii,2},2);
            %set link to 0 for starting
            possibilities{ii,4} = 0;
        end
    end
    %Find single possibilities if they exist
    %First, sort possibilities cell on the size of the possibilities
%     possibilities = sortrows(possibilities,[3 1]);

    %now, have a while loop that continues untill all of the spots have only 1
    %link left
    stuck = 0;
    while (stuck == 0) %% sum([possibilities{:,3}]) > size(possibilities,1) && 
        %Loop for all the possibilities that are not yet linked
        stuck = 0;
        poss_old = possibilities;
        for ii = 1:size(possibilities,1)
            if ii ~= [possibilities{ii,1}]
                disp('Encountered error, look at lines 104-114')
                keyboard
            end
        end
        for ii = find([possibilities{:,3}] > 0)
            %if there's only 1 option which is not yet assigned
            if (possibilities{ii,3} == 1 && possibilities{ii,4} == 0)
%                 currPartner = find([possibilities{:,1}] == [possibilities{ii,2}]);
                currPartner = possibilities{ii,2};
                
                %Set the linking of this and the other to each other
                %Of the one being checked
                possibilities{ii,4} = possibilities{ii,2};
                
                %Set the possibility of the other to this option as well
                %and of its partner
                possibilities{currPartner,2} = possibilities{ii,1};
                possibilities{currPartner,4} = possibilities{ii,1};

                %remove the same option from others
                %For this, check the other arrays
                for k = 1:size(possibilities,1)
                    %check if k is not the current instance
                    if k ~= ii
                        %or its partner
                        if k ~= currPartner
                            %skip this if its empty
                            if ~isempty([possibilities{k,2}])
%                                 check if it will be found
                                if (sum(possibilities{k,2}==ii) + sum(possibilities{k,2}==possibilities{ii,2}))>0

                                    %now, only all others are found
                                    %remove index of i
                                    possibilities{k,2}(find([possibilities{k,2}] == possibilities{ii,1})) = [];
                                    %remove index of partner
                                    possibilities{k,2}(find([possibilities{k,2}] == possibilities{ii,2})) = [];
                                end
                            end
                        end
                    end
                end

                %Sort possibilities again
                for k = 1:size(possibilities,1)
                    possibilities{k,3} = size(possibilities{k,2},2);
                end
                ii = size(possibilities,1)+1; %end the parent loop, go back to the while-loop
            end
        end
        %If nothing changes in a loop, just assume that all other particles are
        %unresolvable.
        if isequal(possibilities, poss_old)
            %if we don't want to guess anything
            if guess_possible_linking == 0
                %Collapse all multiple option possibilities to no option at all
                %(ignored basically)
                for k = 1:size(possibilities,1)
                    if size(possibilities{k,2},2)>1
                        possibilities{k,2} = [];
                    end
                end
                %Get out of the while-loop
                stuck = 1;
            %otherwise, if we want to choose the option that's closest to
            %centre of guesstimate or something...
            else
                %todo
            end
        end
    end
    
    %% Create final x,y positions
    %First, find only 1 original id of linked particles
    %For this, check if the id of the linked particle is smaller than the
    %id of its linked particle. Only keep those ids.
    idlist = [];
    %Loop over all input particles
    for k = 1:size(possibilities,1)
        %check if it's linked
        if ~isempty(possibilities{k,2})
            %check if it's the smallest
            if possibilities{k,1} < possibilities{k,2}
                idlist = [idlist; [k possibilities{k,2}]];
            end
        end
    end
    %Now, create the particles table with x,y,rotation,distance
    particlesT = zeros(size(idlist,1),12);
    for k = 1:size(idlist,1)
        %xypos of linked particles
        xpos1 = loc_array(idlist(k,1),2);
        xpos2 = loc_array(idlist(k,2),2);
        ypos1 = loc_array(idlist(k,1),3);
        ypos2 = loc_array(idlist(k,2),3);
        %other calculations
        particlesT(k,1) = frame;
        particlesT(k,2) = xpos1; %x pos
        particlesT(k,3) = ypos1; %y pos
        particlesT(k,4) = all_array(idlist(k,1),4); %z pos of spatial
        particlesT(k,5) = atan2((ypos2-ypos1),(xpos2-xpos1)); %angle = atan2(dy,dx)
        particlesT(k,6) = sqrt((xpos2-xpos1)^2+(ypos2-ypos1)^2); %distance
        particlesT(k,7) = all_array(idlist(k,1),8); %sigma1 of spatial
        particlesT(k,8) = all_array(idlist(k,1),9); %sigma2 of spatial
        particlesT(k,9) = all_array(idlist(k,2),8); %sigma1 of spectral
        particlesT(k,10) = all_array(idlist(k,2),9); %sigma2 of spectral
        particlesT(k,11) = all_array(idlist(k,1),6); %locuncertainty of spatial
        particlesT(k,12) = all_array(idlist(k,1),7); %intensity of spatial
    end
    particles(curparticle_id:curparticle_id-1+size(particlesT,1),:) = particlesT;
    curparticle_id = curparticle_id+size(particlesT,1);
    
    allpossibilities{frame,:,:} = possibilities;
end
particles(particles(:,1)==0,:) = [];

%% end function
end