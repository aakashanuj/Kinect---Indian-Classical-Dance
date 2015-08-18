%folders = {'N:\Final_MTP_data\utsanga\1\Dancer1\','N:\Final_MTP_data\utsanga\1\Dancer2\','N:\Final_MTP_data\utsanga\1\Dancer3\','N:\Final_MTP_data\tirmana\1\Dancer1\','N:\Final_MTP_data\tirmana\2\Dancer1\','N:\Final_MTP_data\tei_tei_dhatta\1\Dancer1\','N:\Final_MTP_data\tei_tei_dhatta\2\Dancer1\', 'N:\Final_MTP_data\tei_tei_dhatta\3\Dancer1\', 'N:\Final_MTP_data\sarika\1\Dancer1\','N:\Final_MTP_data\sarika\3\Dancer1\','N:\Final_MTP_data\sarika\4\Dancer1\','N:\Final_MTP_data\pakka\1\Dancer1\','N:\Final_MTP_data\pakka\2\Dancer1\','N:\Final_MTP_data\pakka\3\Dancer1\','N:\Final_MTP_data\pakka\4\Dancer1\','N:\Final_MTP_data\paikkal\1\Dancer1\','N:\Final_MTP_data\paikkal\3\Dancer1\','N:\Final_MTP_data\natta\1\Dancer1\'}
folders = {'P:\Final_MTP_data\utsanga\1\Dancer1\', 'P:\Final_MTP_data\tirmana\1\Dancer1\', 'P:\Final_MTP_data\tei_tei_dhatta\1\Dancer1\', 'P:\Final_MTP_data\sarika\1\Dancer1\', 'P:\Final_MTP_data\pakka\1\Dancer1\', 'P:\Final_MTP_data\paikkal\1\Dancer1\', 'P:\Final_MTP_data\natta\1\Dancer1\', 'P:\Final_MTP_data\joining\1\Dancer2\', 'P:\Final_MTP_data\katti_kartari\1\Dancer2\', 'P:\Final_MTP_data\kuditta_nattal\1\Dancer2\', 'P:\Final_MTP_data\tatta\3\Dancer1\', 'P:\Final_MTP_data\kuditta_mettu\3\Dancer1\', 'P:\Final_MTP_data\kuditta_tattal\1\Dancer2\', 'P:\Final_MTP_data\sarrikkal\1\Dancer1\','P:\Final_MTP_data\mandi\1\Dancer1\'}
outFileName = strcat('featureVectors.txt');
fid = fopen(outFileName,'wt');
framesCounter = 1;
framesToCluster{framesCounter} = 'temporary'; 
startingIndex = 0;
indicesOfAdavus = zeros(1,30);
indexCount = 1;
for folderCount = 1:size(folders, 2)
%for folderCount = 1:8
    featureCountForThisAdavu = 0;
    audiofiles = dir(strcat(folders{folderCount},'*.wav'));
    o = mironsets(strcat(folders{folderCount},audiofiles(1).name));
    
    % Removing peaks which are very close
    times = mirgetdata(o)
    times = times';
    d = get(o,'PeakVal')
    values = d{1}{1}{1}
    a = zeros(1,size(times,2))
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
            clearvars -except indicesOfAdavus indexCount startingIndex endIndex featureCountForThisAdavu framesToCluster framesCounter array folders folderCount j e fid;
            continue;
        end

        framesToCluster{framesCounter} = strcat(folders{folderCount},'color_USB-VID_045E&PID_02BF-0000000000000000_',num2str(j),'.png');
        framesCounter = framesCounter + 1;
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
        %featureAttributeCount = 1;
        %featureTemp = zeros(1,34);
        for i = 1:count-2
            % Ignore the hip center joint, since it has NaN and NaN
            if (i==1 || i==16)
                continue;
            end
            if i==19
                fprintf(fid,'%f\t%f\n',thetaY{i},thetaXZ{i});
                %featureTemp(featureAttributeCount) = thetaY{i};
                %featureTemp(featureAttributeCount + 1) = thetaXZ{i};
                %featureAttributeCount = featureAttributeCount + 2;
            else
                fprintf(fid,'%f\t%f\t',thetaY{i},thetaXZ{i});
                %featureTemp(featureAttributeCount) = thetaY{i};
                %featureTemp(featureAttributeCount + 1) = thetaXZ{i};
                %featureAttributeCount = featureAttributeCount + 2;
            end
        end
        clearvars -except indicesOfAdavus indexCount startingIndex endIndex featureCountForThisAdavu framesToCluster framesCounter  folders array e folderCount j fid;
    end
    disp(featureCountForThisAdavu);
    indicesOfAdavus(indexCount) = startingIndex + 1;
    indicesOfAdavus(indexCount+1) = startingIndex + featureCountForThisAdavu;
    startingIndex = indicesOfAdavus(indexCount+1);
    indexCount = indexCount+ 2;
    clearvars -except indicesOfAdavus indexCount startingIndex endIndex framesToCluster framesCounter folders fid folderCount;
    %figure,imshow(strcat('G:\backup\',folders{folderCount},'\color_USB-VID_045E&PID_02BF-0000000000000000_',num2str(best_posture),'.png'))
end

fclose(fid);
fid = fopen('featureVectors.txt');
tline = fgets(fid);
count = 1;
X = zeros(1,34);
while ischar(tline)
    disp(tline);
    parts = regexp(tline,'\t','split');
    for i = 1:size(parts,2)
        %disp(str2double(parts(i)));
        tempArray(i) = str2double(parts(i));
    end
    X(count, :) = tempArray;
    tline = fgets(fid);
    count = count + 1;
end
disp(X);
disp(size(X));

save('firstStepData')
disp('Saving done');
pause;