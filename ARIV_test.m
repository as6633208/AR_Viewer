function varargout = ARIV_test(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 16-Jun-2017 14:06:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ARIV_test_OpeningFcn, ...
                   'gui_OutputFcn',  @ARIV_test_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function ARIV_test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axis off;

function varargout = ARIV_test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Camera 1 On.
function pushbutton2_Callback(hObject, eventdata, handles)
vid1 = videoinput('winvideo',1);
hImage=image(zeros(640,480,1),'Parent',handles.axes2);
preview(vid1,hImage);

% --- Stop Camera.
function pushbutton10_Callback(hObject, eventdata, handles)
closepreview;

% --- Camera 2 On.
function pushbutton4_Callback(hObject, eventdata, handles)
vid2 = videoinput('winvideo',2);
hImage=image(zeros(640,480,1),'Parent',handles.axes2);
preview(vid2,hImage);

% --- Calibrator.
function pushbutton5_Callback(hObject, eventdata, handles)
% Define images to process

imageFileNames1 = {'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_1.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_10.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_11.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_12.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_15.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_2.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_3.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_4.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_6.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_8.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\L\imgFromCamL_9.png',...
};
imageFileNames2 = {'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_1.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_10.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_11.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_12.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_15.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_2.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_3.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_4.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_6.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_8.png',...
    'C:\opencv\AR\Stereo camera\Fundamental matrix\R\imgFromCamR_9.png',...
};

% Detect checkerboards in images
[imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints(imageFileNames1, imageFileNames2);

% Generate world coordinates of the checkerboard keypoints
squareSize = 24;  % in units of 'mm'
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

% Calibrate the camera
[stereoParams, pairsUsed, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints, ...
'EstimateSkew', false, 'EstimateTangentialDistortion', false, ...
'NumRadialDistortionCoefficients', 2, 'WorldUnits', 'mm', ...
'InitialIntrinsicMatrix', [], 'InitialRadialDistortion', []);

% View reprojection errors
h1=figure; showReprojectionErrors(stereoParams, 'BarGraph');

% Visualize pattern locations
h2=figure; showExtrinsics(stereoParams, 'CameraCentric');

% --- Position.

function pushbutton6_Callback(hObject, eventdata, handles)

clc;
cam1 = videoinput('winvideo', 1, 'RGB24_640x480'); % 1 = 1st camera
cam2 = videoinput('winvideo', 2, 'RGB24_640x480'); % 2 = 2nd camera

% Capture one frame to get its size.
img1 = getsnapshot(cam1); frameSize1 = size(img1);
img2 = getsnapshot(cam2); frameSize2 = size(img2);

load('stereoParams.mat');

[imgRect1, imgRect2] = rectifyStereoImages(img1,img2,stereoParams);

%figure;
%imshow(stereoAnaglyph(imgRect1, imgRect2)); title('Rectified Image');

imgGray1 = rgb2gray(imgRect1);
imgGray2 = rgb2gray(imgRect2);

disparityMap = disparity(imgGray1,imgGray2);
%figure;
%imshow(disparityMap, [0, 64]); title('Disparity Map');
colormap jet
colorbar

point3D = reconstructScene(disparityMap, stereoParams);

% Limit the range of Z and X for display.
z = point3D(:,:,3);
z(z < 0 | z > 8000) = NaN;
x = point3D(:,:,1);
x(x < -3000 | x > 3000) = NaN;
point3D(:,:,3) = z;
point3D(:,:,1) = x;

%figure;
%pcshow(point3D, imgRect1, 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
set(gca, 'CameraViewAngle',10, 'CameraUpVector',[0 -1 0],...
    'CameraPosition',[16500 -13852 -49597], 'DataAspectRatio',[1 1 1]);
title('Reconstructed 3-D Scene');
numImages = 20;
images1 = cell(1, numImages);
images2 = cell(1, numImages);
for i = 1:numImages
    images1{i} = fullfile('c:\','opencv','AR','Stereo camera','Fundamental matrix','L',sprintf('imgFromCamL_%d.png',i));
    images2{i} = fullfile('c:\','opencv','AR','Stereo camera','Fundamental matrix','R',sprintf('imgFromCamR_%d.png',i));
end

cameraParams = cameraParameters('IntrinsicMatrix',stereoParams.CameraParameters1.IntrinsicMatrix,'RadialDistortion',stereoParams.CameraParameters1.RadialDistortion);
J = undistortImage(img1,cameraParams);
global LOC_OBJ_ NUM_OF_OBJ DIS_OF_OBJ_;
NUM_OF_OBJ=1;
imshow(J,'Parent',handles.axes2);
while NUM_OF_OBJ>=0
    objectRegion = round(getPosition(imrect)); % mouse click and drag
    confirm = questdlg('Is the area you choose correct?','Confirm','Yes','No','');
    switch confirm
        case 'Yes'
            if NUM_OF_OBJ == 1
                % Show object region with a bounding box in the first camera.
                object = insertShape(J, 'Rectangle', objectRegion, 'LineWidth',2,'Color','yellow');
                object2=insertText(object,objectRegion(1:2),NUM_OF_OBJ,'FontSize',12,'BoxColor','red','TextColor','white');
                imshow(object2,'Parent',handles.axes2); title('Red Box Shows Object Region');
                
                   % Find the centroids of object. 
                centroids = [round(objectRegion(:, 1) + objectRegion(:,3 ) / 2),... % we want a digital image pixel value, so we round to the nearest row and column
                round(objectRegion(:, 2) + objectRegion(:, 4) / 2)];

                % Find the 3-D world coordinates of the centroids.
                centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
                X = point3D(:, :, 1);
                Y = point3D(:, :, 2);
                Z = point3D(:, :, 3);
                centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];

                % Display detected object with 3-D world coordinates.
                reshapeCentroids3D = reshape(centroids3D,[1,3]); % reshape into 1 row, 3 column
               %object2 = insertObjectAnnotation(object2,'rectangle', objectRegion,...
                %['XYZ: ',num2str(reshapeCentroids3D)], 'FontSize', 18);
            
                LOC_OBJ_{NUM_OF_OBJ}=num2str(reshapeCentroids3D);
                
                %imshow(object2,'Parent',handles.axes2);
                %title('Detected object with 3-D world coordinates (scale in millimeter)');

                % Find the distances from the camera in meters.
                dists = sqrt(sum(centroids3D .^ 2)) / 1000;

                % Display the detected object and distances.
                labels = cell(1, numel(dists));
                for i = 1:numel(dists)
                    labels{i} = sprintf('Distance: %0.2f meters', dists(i));
                end
                DIS_OF_OBJ_{NUM_OF_OBJ}=labels;
                %imshow(insertObjectAnnotation(object, 'rectangle', objectRegion, labels,'FontSize', 18),'Parent',handles.axes2);
                %title('Detected Object with distance (scale in meter)');
                
            else
                % Show object region with a bounding box in the first camera.
                object2 = insertShape(object2, 'Rectangle', objectRegion, 'LineWidth',2,'Color','yellow');
                object2=insertText(object2,objectRegion(1:2),NUM_OF_OBJ,'FontSize',12,'BoxColor','red','TextColor','white');
                imshow(object2,'Parent',handles.axes2); title('Red Box Shows Object Region');
                
                    % Find the centroids of object. 
                centroids = [round(objectRegion(:, 1) + objectRegion(:,3 ) / 2),... % we want a digital image pixel value, so we round to the nearest row and column
                round(objectRegion(:, 2) + objectRegion(:, 4) / 2)];

                % Find the 3-D world coordinates of the centroids.
                centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
                X = point3D(:, :, 1);
                Y = point3D(:, :, 2);
                Z = point3D(:, :, 3);
                centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];

                % Display detected object with 3-D world coordinates.
                reshapeCentroids3D = reshape(centroids3D,[1,3]); % reshape into 1 row, 3 column
                %object2 = insertObjectAnnotation(object2,'rectangle', objectRegion,...
                %['XYZ: ',num2str(reshapeCentroids3D)], 'FontSize', 18);
                
                LOC_OBJ_{NUM_OF_OBJ}=num2str(reshapeCentroids3D);
            
                %imshow(object2,'Parent',handles.axes2);
                %title('Detected object with 3-D world coordinates (scale in millimeter)');

                % Find the distances from the camera in meters.
                dists = sqrt(sum(centroids3D .^ 2)) / 1000;

                % Display the detected object and distances.
                labels = cell(1, numel(dists));
                for i = 1:numel(dists)
                    labels{i} = sprintf('Distance: %0.2f meters', dists(i));
                end
                DIS_OF_OBJ_{NUM_OF_OBJ}=labels;
                %imshow(insertObjectAnnotation(object, 'rectangle', objectRegion, labels,'FontSize', 18),'Parent',handles.axes2);
                %title('Detected Object with distance (scale in meter)');
                
            end
            NUM_OF_OBJ=NUM_OF_OBJ+1;
 
        case 'No'
            if NUM_OF_OBJ == 1
                NUM_OF_OBJ=1;
            end
    end
end


% --- Choose Data
function pushbutton7_Callback(hObject, eventdata, handles)
global LOC_OBJ_ NUM_OF_OBJ DIS_OF_OBJ_ DATA_OF_OBJ_;
set(handles.text4,'visible','on')
set(handles.edit3,'visible','on')
set(handles.pushbutton11,'visible','on')
set(handles.popupmenu1,'visible','on')

for i = 1:NUM_OF_OBJ-1
    list{i}=i;
end
set(handles.popupmenu1,'String',list);


% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global LOC_OBJ_ NUM_OF_OBJ DIS_OF_OBJ_ DATA_OF_OBJ_;
set(handles.text4,'visible','off')
set(handles.edit3,'visible','off')
set(handles.pushbutton11,'visible','off')
set(handles.popupmenu1,'visible','off')
set(handles.text5,'visible','on')
set(handles.text7,'visible','on')

for i = 1:NUM_OF_OBJ-1
    list_of_data{i}=['Data of ',int2str(i),' is ',DATA_OF_OBJ_{i},'.    Location of ',int2str(i),' is ',LOC_OBJ_{i},'.    Distant of ',int2str(i),' is ',DIS_OF_OBJ_{i}];
    
end
set(handles.text5,'String',list_of_data);






% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
%==========================================================================
% cam1
%==========================================================================
cam1 = videoinput('winvideo', 1, 'RGB24_640x480');
set(cam1, 'ReturnedColorSpace', 'RGB');
imgFromCam1 = getsnapshot(cam1);
 
% This is where and what my image will be saved
counter  = 1;
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\L\';
baseName = 'imgFromCamL_';
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
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\R\';
baseName = 'imgFromCamR_';
newName  = [baseDir baseName num2str(counter) '.png'];
while exist(newName,'file')
    counter = counter + 1;
    newName = [baseDir baseName num2str(counter) '.png'];
end    
imwrite(imgFromCam2, newName);
 
% Display image from cam1 and cam2.
imshowpair(imgFromCam1, imgFromCam2,'montage', 'Parent', handles.axes2);
% Clear
clear all
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DATA_OF_OBJ_
value=get(handles.popupmenu1,'Value');
DATA_OF_OBJ_{value}=get(handles.edit3,'String');
set(handles.edit3,'String','');



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
cam = ipcam('http://192.168.0.198:80/mjpg/video.mjpg','root','0000');
hImage=image(zeros(640,480,1),'Parent',handles.axes2);
preview(cam);
imshow(snapshot(cam),'Parent',handles.axes2);


