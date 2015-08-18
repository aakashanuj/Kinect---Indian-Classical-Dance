% This contains the list of folders which have the dance data extracted in
% them
folders = {'P:\Final_MTP_data\tirmana\1\Dancer1\','P:\Final_MTP_data\utsanga\1\Dancer3\'};
for folderCount = 1 : size(folders,2)
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
        endUUt
    end
    e = goodtimes;
    e = e*30;
    e = int32(e);

    outputVideo = VideoWriter(strcat('aakash',num2str(folderCount),'.flv'));
    outputVideo.FrameRate = 2;
    open(outputVideo)
    for i = 1 : size(e,2)                       
        img = imread(strcat(folders{folderCount},'color_USB-VID_045E&PID_02BF-0000000000000000_',num2str(e(i)),'.png'));
        writeVideo(outputVideo,img);
    end
    close(outputVideo)
end
