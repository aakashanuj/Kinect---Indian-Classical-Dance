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
imagen = frame.cdata;

    % Save the image to disk
    imwrite(imagen, 'aa.png');
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

hf=figure; 

%// 1st subplot
hs1 = subplot(2,1,1);
imshow(a);

%// Create uipanel and place it below current image/subplot 1

hPanel= uipanel('Units','Normalized');

% Set the uipanel position to where the image is desired to be displayed.

set(hPanel,'position',[0 0.05 1 .515]);

%// Create an axes whose parent is the uipanel.
ax1 = axes('parent',hPanel,'position',[0 0 1 1],'Units','normalized');

%// Display image to get the handle.
him1 = imshow(im,'parent',ax1);

% Create the scroll panel
hSP = imscrollpanel(hPanel,him1);
api = iptgetapi(hSP);
api.setVisibleLocation(640,0.5);