%==========================================================================
% cam1
%==========================================================================
cam1 = videoinput('winvideo', 1, 'RGB24_640x480');
set(cam1, 'ReturnedColorSpace', 'RGB');
imgFromCam1 = getsnapshot(cam1);
 
% This is where and what my image will be saved
counter  = 1;
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\L\';
baseName = 'imgFromCam1_';
newName  = [baseDir baseName num2str(counter) '.png'];
while exist(newName,'file')
    counter = counter + 1;
    newName = [baseDir baseName num2str(counter) '.png'];
end    
imwrite(imgFromCam1, newName);
 
%==========================================================================
% cam2
%==========================================================================
cam2 = videoinput('winvideo', 2, 'RGB24_640x480');
set(cam2, 'ReturnedColorSpace', 'RGB');
imgFromCam2 = getsnapshot(cam2);
 
% This is where and what my image will be saved
counter  = 1;
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\L\';
baseName = 'imgFromCam2_';
newName  = [baseDir baseName num2str(counter) '.png'];
while exist(newName,'file')
    counter = counter + 1;
    newName = [baseDir baseName num2str(counter) '.png'];
end    
imwrite(imgFromCam2, newName);
 
% Display image from cam1 and cam2.
figure
imshowpair(imgFromCam1, imgFromCam2,'montage');
 
% Clear
clear all
