clear
clc

% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
imagefiles = dir('D:\temp\color*.png');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(strcat('D:\temp\',currentfilename));
   images{ii} = currentimage;
end
disp(nfiles);

audiofiles = dir('D:\temp\*.wav');      
o = mironsets(strcat('D:\temp\',audiofiles(1).name))
frame = getframe(gcf);
dcm_obj = datacursormode(gcf);
set(dcm_obj,'DisplayStyle','datatip',...
'SnapToDataVertex','off','Enable','on')

d = get(o,'PeakPos')
e = d{1}{1}{1}
e = e /100.0
e = e*30
e = int32(e)

a = imread('aa.png');
textColor    = [255, 255, 255]; % [red, green, blue]
textLocation = [240 50];       % [x y] coordinates
textInserter = vision.TextInserter('Selected !!', 'Color', textColor, 'FontSize', 40, 'Location', textLocation);
for ii=1:nfiles 
    flag =0;
    for j = 1:size(e, 2)
        if e(j)==ii
            flag =1;
        end
    end
    if flag==1
        images{ii} = step(textInserter, images{ii});
    end
end
im = cat(2,images{:});

hFig = figure('Toolbar','none',...
              'Menubar','none');
hIm = imshow(im);
hSP = imscrollpanel(hFig,hIm);
set(hSP,'Units','normalized',...
        'Position',[0 .1 1 .9])
api = iptgetapi(hSP);
api.setVisibleLocation(640,0.5);

while(1)
    pause;
    c_info = getCursorInfo(dcm_obj);
    xval  = c_info.Position(1);
    api.setVisibleLocation(640 * 30 * xval - (640*3),0.5);
end