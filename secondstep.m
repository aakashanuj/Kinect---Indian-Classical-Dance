% X contains the feature vectors
% Now we try clustering
load('firstStepData')

DOKMEANS = 0;
if DOKMEANS == 1
    INITIAL_NUMBER_OF_CLUSTERS = 10;
    FINAL_NUMBER_OF_CLUSTERS = 100;
    NUMBER_OF_SIMULATIONS = 1;
    Clusters = zeros(1, NUMBER_OF_SIMULATIONS);
    means = zeros(1,FINAL_NUMBER_OF_CLUSTERS - INITIAL_NUMBER_OF_CLUSTERS + 1);
    x = INITIAL_NUMBER_OF_CLUSTERS:FINAL_NUMBER_OF_CLUSTERS;
    for j = 1:NUMBER_OF_SIMULATIONS % Number of simulations
        optimalNumberOfClusters = -1;
        maxMean = -1;
        for i = INITIAL_NUMBER_OF_CLUSTERS:FINAL_NUMBER_OF_CLUSTERS
            disp(i);
            %idx3 =  KMeansCustom(X,i);
            idx3 =  kmeans(X,i);
            %[silh3,h] = silhouette(X,idx3,@distanceFunction);
            [silh3,h] = silhouette(X,idx3);
            meanValue = mean(silh3);
            %meanValue = size(find(silh3>0.4),1);
            means(i-(INITIAL_NUMBER_OF_CLUSTERS-1)) = meanValue;
            disp(meanValue);
            if meanValue > maxMean
                maxMean = meanValue;
                optimalNumberOfClusters = i;
            end
        end
        disp(optimalNumberOfClusters);
        Clusters(j) = optimalNumberOfClusters;
    end
    figure, plot(x, means)
end
%optimalNumberOfClusters = mode(Clusters)
optimalNumberOfClusters = 74;
pause;

% Now we know the optimal number of clusters for K means, go ahead and do K means now
C = zeros(size(X,1),34);
Z = zeros(size(X,1),1);
for j = 1:1 % Number of simulations
    %idx3 =  KMeansCustom(X,optimalNumberOfClusters);
    [idx3, C] =  kmeans(X,optimalNumberOfClusters);
    Z(:,j) = idx3;
end

% Now we want to assign a class to each of the frames extracted
for y = 1: size(Z,1)
    Y(y) = mode(Z(y,:));
end
disp(Y)

%pause;
disp('About to write images');

WRITE = 0;
if(WRITE==1)
    %Now Y contains the class for each data point
    %Create a separate image directory for each cluster
    for cluster = 1:optimalNumberOfClusters
        folderName = strcat('Cluster',num2str(cluster));
        mkdir(folderName);
        for k = 1:size(Z,1)
            if Y(k) == cluster
                tempstr = cellstr(framesToCluster{k});
                image = imread(strrep(tempstr{1},'N:','P:'));
                imwrite(image,strcat('G:\backup\',folderName,'\',num2str(k),'.png'),'png');
            end
        end
    end
end

% The list of Adavus for which we have proper sound files
adavus{1} = 'Utsanga';
adavus{2} = 'Tirmana';
adavus{3} = 'Tei Tei Dhatta';
adavus{4} = 'Sarika';
adavus{5} = 'Pakka';
adavus{6} = 'Paikkal';
adavus{7} = 'Natta';
adavus{8} = 'Joining';
adavus{9} = 'Katti/Kartari';
adavus{10} = 'Kuditta Nattal';
adavus{11} = 'Tatta';
adavus{12} = 'Kuditta Mettu';
adavus{13} = 'Kuditta Tattal';
adavus{14} = 'Sarrikkal';
adavus{15} = 'Mandi';

vecCount = 1;
for i = 1:2:size(indicesOfAdavus,2)
    vectorOfCodes{vecCount} = idx3(indicesOfAdavus(i):indicesOfAdavus(i+1));
    vecCount = vecCount + 1;
end
disp(vectorOfCodes);
% Now the vectorOfCodes contains the posture labels for the key postures in
% every adavu
% We should perform edit distance algo on them with key postures of the new adavu


featuresNew{1} = [0 0 0];
%folders = {'N:\Final_MTP_data\utsanga\1\Dancer1\', 'N:\Final_MTP_data\tirmana\1\Dancer1\', 'N:\Final_MTP_data\tei_tei_dhatta\1\Dancer1\', 'N:\Final_MTP_data\sarika\1\Dancer1\', 'N:\Final_MTP_data\pakka\1\Dancer1\', 'N:\Final_MTP_data\paikkal\1\Dancer1\', 'N:\Final_MTP_data\natta\1\Dancer1\', 'N:\Final_MTP_data\joining\1\Dancer2\', 'N:\Final_MTP_data\katti_kartari\1\Dancer2\', 'N:\Final_MTP_data\kuditta_nattal\1\Dancer2\'}
folders = {'P:\Final_MTP_data\utsanga\1\Dancer2\', 'P:\Final_MTP_data\tirmana\1\Dancer2\', 'P:\Final_MTP_data\tei_tei_dhatta\1\Dancer3\','P:\Final_MTP_data\sarika\1\Dancer2\', 'P:\Final_MTP_data\pakka\1\Dancer2\','P:\Final_MTP_data\paikkal\1\Dancer2\', 'P:\Final_MTP_data\natta\1\Dancer2\','P:\Final_MTP_data\joining\1\Dancer3\', 'P:\Final_MTP_data\katti_kartari\1\Dancer3\', 'P:\Final_MTP_data\kuditta_nattal\1\Dancer3\','P:\Final_MTP_data\tatta\3\Dancer2\', 'P:\Final_MTP_data\kuditta_mettu\3\Dancer2\', 'P:\Final_MTP_data\kuditta_tattal\1\Dancer3\', 'P:\Final_MTP_data\sarrikkal\1\Dancer1\','P:\Final_MTP_data\mandi\1\Dancer1\'};
%for folderCount = 1:size(folders,2)
for folderCount = 1 : size(folders,2)
    featureCountForThisAdavu = 0;
    audiofiles = dir(strcat(folders{folderCount},'*.wav'));
    o = mironsets(strcat(folders{folderCount},audiofiles(1).name));
    
    % Removing peaks which are very close
    times = mirgetdata(o);
    times = times';
    d = get(o,'PeakVal')
    values = d{1}{1}{1};
    a = zeros(1,size(times,2));
    for i = 1:size(times,2)
        for j = i+1:size(times,2)
            if(abs(times(i) - times(j))<0.15)
                if(values(i)>values(j))
                    a(j)=1;
                end
                if(values(j)>values(i))
                    a(i)=1;
                end
            end
        end
    end

    peakcount = 1;
    goodtimes = zeros(1, size(a,2) - sum(a));
    goodpeaks = zeros(1, size(a,2) - sum(a));
    for i = 1:size(a,2)
        if a(i)==0
            goodtimes(peakcount) = times(i);
            goodpeaks(peakcount) = values(i);
            peakcount = peakcount +1;
        end
    end
    e = goodtimes;
    e = e*30;
    e = int32(e);

    disp(folders{folderCount});

    dirName = strcat(folders{folderCount},'*.mat');
    allFiles = dir(dirName);
    filenames = {allFiles(~[allFiles.isdir]).name};

    for j = 1:size(filenames,2)
        % Check if it is a required frame
        flag =0;
        for k = 1:size(e, 2)
            if e(k)==j
                flag =1;
            end
        end
        if flag==0
            continue;
        end
        %Checking done

        matFileToLoad = strcat(folders{folderCount},'USB-VID_045E&PID_02BF-0000000000000000_',num2str(j),'.mat');
        load(matFileToLoad);

        % Find the skeleton which is tracked
        skeleton = -1;
        skeletonAlreadyFound = 0;
        for i = 1:6
            if (strcmp(SkeletonFrame.Skeletons(i).TrackingState,'Tracked')~=0)
                if(skeleton ~= -1)
                    skeletonAlreadyFound = 1;
                end
                skeleton = i;
            end
        end

        % Handle the case if skeleton is not tracked
        if (skeleton == -1 || skeletonAlreadyFound==1)
            for i = 1:19
                % Ignore the hip center joint, since it has NaN and NaN
                if (i==1 || i==16)
                    continue;
                end
                if i==19
                    %fprintf(fid,'%f\t%f\n',1000000,1000000);
                    fault = 1;
                else
                    %fprintf(fid,'%f\t%f\t',1000000,1000000);
                    fault = 1;
                end
            end
            clearvars -except X C idx3 vectorOfCodes featuresNew adavus indicesOfAdavus indexCount startingIndex endIndex featureCountForThisAdavu framesToCluster framesCounter array folders folderCount j e fid;
            continue;
        end

        % Data structures to store the outputs
        jointNames={};
        thetaY = {};
        thetaXZ ={};
        count = 1;

        % Translate to origin
        bx10=SkeletonFrame.Skeletons(skeleton).Joints(1).Position.X;
        by10=SkeletonFrame.Skeletons(skeleton).Joints(1).Position.Y;
        bz10=SkeletonFrame.Skeletons(skeleton).Joints(1).Position.Z;
        for iter=1:20
            SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.X = SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.X-bx10;
            SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.Y = SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.Y-by10;
            SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.Z = SkeletonFrame.Skeletons(skeleton).Joints(iter).Position.Z-bz10;
        end  

        % Consider all 20 bones of the skeleton found
        for bone = 1:20
            startIndex = SkeletonFrame.Skeletons(skeleton).BoneOrientations(bone).StartJointIndex;
            endIndex = SkeletonFrame.Skeletons(skeleton).BoneOrientations(bone).EndJointIndex;
            startPos = SkeletonFrame.Skeletons(skeleton).Joints(startIndex).Position;
            endPos = SkeletonFrame.Skeletons(skeleton).Joints(endIndex).Position;
            parentVector = [ SkeletonFrame.Skeletons(skeleton).Joints(startIndex).Position.X SkeletonFrame.Skeletons(skeleton).Joints(startIndex).Position.Y SkeletonFrame.Skeletons(skeleton).Joints(startIndex).Position.Z];
            childVector = [ SkeletonFrame.Skeletons(skeleton).Joints(endIndex).Position.X SkeletonFrame.Skeletons(skeleton).Joints(endIndex).Position.Y SkeletonFrame.Skeletons(skeleton).Joints(endIndex).Position.Z];
            % For each bone, subtract the parent joint position from the child joint position so that a vector is formed that crosses through the origin of the coordinate system. Normalize this vector, lets call it diff.
            diffVector_orig = childVector - parentVector;
            diffVector = diffVector_orig / norm(diffVector_orig);
            matrix = SkeletonFrame.Skeletons(skeleton).BoneOrientations(1).AbsoluteRotation.Matrix(1:3,1:3);
            diffVector = inv(matrix) * diffVector';

            % Calculating the angle with the Y axis, with orientation independence
            y =  [0, 1, 0];
            dotProduct = dot(y,diffVector);
            angle_y = acosd(dotProduct);
            thetaY{count} = angle_y;

            % Calculating the angle with the x-z plane
            x = [1,0];
            vector2D = [diffVector(1) diffVector(3)];
            vector2D = vector2D / norm(vector2D);
            dotProduct = dot(x, vector2D);
            angle_xz = acosd(dotProduct);
            thetaXZ{count} = angle_xz;

            jointNames{count}=strcat(SkeletonFrame.Skeletons(skeleton).Joints(startIndex).JointType,SkeletonFrame.Skeletons(skeleton).Joints(endIndex).JointType);
            count = count + 1;
            %disp(strcat(SkeletonFrame.Skeletons(skeleton).Joints(startIndex).JointType,SkeletonFrame.Skeletons(skeleton).Joints(endIndex).JointType));
        end

        featureCountForThisAdavu = featureCountForThisAdavu + 1;
        featureAttributeCount = 1;
        featureTemp = zeros(1,34);
        for i = 1:count-2
            % Ignore the hip center joint, since it has NaN and NaN
            if (i==1 || i==16)
                continue;
            end
            if i==19
                %fprintf(fid,'%f\t%f\n',thetaY{i},thetaXZ{i});
                featureTemp(featureAttributeCount) = thetaY{i};
                featureTemp(featureAttributeCount + 1) = thetaXZ{i};
                featureAttributeCount = featureAttributeCount + 2;
            else
                %fprintf(fid,'%f\t%f\t',thetaY{i},thetaXZ{i});
                featureTemp(featureAttributeCount) = thetaY{i};
                featureTemp(featureAttributeCount + 1) = thetaXZ{i};
                featureAttributeCount = featureAttributeCount + 2;
            end
        end
        clearvars -except X C idx3 featuresNew featureTemp vectorOfCodes adavus  indicesOfAdavus indexCount startingIndex endIndex featureCountForThisAdavu framesToCluster framesCounter  folders array e folderCount j fid;
        % Now find the closest cluster for this feature vector
        featuresNew{featureCountForThisAdavu} = featureTemp;
    end
    
    %Check 
    B = zeros(featureCountForThisAdavu,34);
    disp(size(B));
    for iter = 1 :featureCountForThisAdavu
        B(iter,:) = featuresNew{iter};
    end
    
    classesNew = knnclassify(B,X,idx3);
    
    index_min_standard = -1;
    index_min_weighted = -1;
    minStandardEditDistance = 100000;
    minWeightedEditDistance = 100000;
    
    index_min_standard1 = -1;
    index_min_weighted1 = -1;
    minStandardEditDistance1 = 100000;
    minWeightedEditDistance1 = 100000;
    
    index_max_standard = -1;
    index_max_weighted = -1;
    maxStandardEditDistance = -100000;
    maxWeightedEditDistance = -100000;
    
    index_max_standard1 = -1;
    index_max_weighted1 = -1;
    maxStandardEditDistance1 = -100000;
    maxWeightedEditDistance1 = -100000;
    
    for iter = 1:size(vectorOfCodes,2)
        disp(adavus{iter});
        
        tempsize = size(vectorOfCodes{iter}',2);
        tempsize1 = size(classesNew',2);
        
        newVec = ones(1,tempsize);
        newVec = newVec * 75;
        
        newVec1 = ones(1,tempsize1);
        newVec1 = newVec1 * 75;
        
        % This vector is vertical, make it horizontal later
        concatVec = cat(1,vectorOfCodes{iter} , newVec');
        concatVec1 = cat(1, classesNew, newVec1');
        
        edit = findEditDistance(concatVec',concatVec1');
        disp(edit);
        if edit<minStandardEditDistance
            minStandardEditDistance = edit;
            index_min_standard = iter;
        elseif edit<minStandardEditDistance1
            minStandardEditDistance1 = edit;
            index_min_standard1 = iter;
        end
        
        if edit>maxStandardEditDistance
            maxStandardEditDistance = edit;
            index_max_standard = iter;
        elseif edit>maxStandardEditDistance1
            maxStandardEditDistance1 = edit;
            index_max_standard1 = iter;
        end
        
        edit1 = findWeightedEditDistance(concatVec',concatVec1', C);
        disp(edit1);
        if edit1<minWeightedEditDistance
            minWeightedEditDistance = edit1;
            index_min_weighted = iter;
        elseif edit1<minWeightedEditDistance1
            minWeightedEditDistance1 = edit1;
            index_min_weighted1 = iter;
        end
        
        if edit1>maxWeightedEditDistance
            maxWeightedEditDistance = edit1;
            index_max_weighted = iter;
        elseif edit1>maxWeightedEditDistance1
            maxWeightedEditDistance1 = edit1;
            index_max_weighted1 = iter;
        end
        
    end
    
    %disp(strcat('Standard Closest to ', adavus{index_min_standard},' ', adavus{index_min_standard1}, ' ', num2str(minStandardEditDistance), ' ', num2str(minStandardEditDistance1)));
    %disp(strcat('Weighted Closest to ', adavus{index_min_weighted},' ', adavus{index_min_weighted1}, ' ', num2str(minWeightedEditDistance), ' ', num2str(minWeightedEditDistance1)));
    
    disp(sprintf('Standard Closest to %s %s %s %s',adavus{index_min_standard},adavus{index_min_standard1}, num2str(minStandardEditDistance), num2str(minStandardEditDistance1)));
    disp(sprintf('Weighted Closest to %s %s %s %s',adavus{index_min_weighted},adavus{index_min_weighted1}, num2str(minWeightedEditDistance), num2str(minWeightedEditDistance1)));
    
    disp(sprintf('Standard farthest to %s %s %s %s',adavus{index_max_standard},adavus{index_max_standard1}, num2str(maxStandardEditDistance), num2str(maxStandardEditDistance1)));
    disp(sprintf('Weighted farthest to %s %s %s %s',adavus{index_max_weighted},adavus{index_max_weighted1}, num2str(maxWeightedEditDistance), num2str(maxWeightedEditDistance1)));
    clearvars -except X C idx3 featuresNew vectorOfCodes adavus indicesOfAdavus indexCount startingIndex endIndex framesToCluster framesCounter folders fid folderCount;

    %figure,imshow(strcat('G:\backup\',folders{folderCount},'\color_USB-VID_045E&PID_02BF-0000000000000000_',num2str(best_posture),'.png'))
end






