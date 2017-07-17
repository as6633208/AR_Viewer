function varargout = untitled(varargin)
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

% Last Modified by GUIDE v2.5 16-Jul-2017 12:08:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
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
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --- 相機右啟動
function pushbutton2_Callback(hObject, eventdata, handles)
%cam_R=videoinput('winvideo',2);
%hImage=image(zeros(800,1000,1),'Parent',handles.axes2);
%preview(cam_R,hImage);
Cam_R = imread('C:\FFOutput\imgFromCamR_1.png');
imshow(Cam_R,'Parent',handles.axes2);
% --- 停止當前相機
function pushbutton10_Callback(hObject, eventdata, handles)
closepreview;
% --- 相機左啟動
function pushbutton4_Callback(hObject, eventdata, handles)
%cam_L=videoinput('winvideo',1);
%hImage=image(zeros(800,1000,1),'Parent',handles.axes2);
%preview(cam_L,hImage);
Cam_L = imread('C:\FFOutput\imgFromCamL_1.png');
imshow(Cam_L,'Parent',handles.axes2);
% --- 截圖
function pushbutton9_Callback(hObject, eventdata, handles)
%==========================================================================
% cam1
%==========================================================================
%cam1 = videoinput('winvideo', 1, 'RGB24_960x720');
%Cam_L = ipcam('http://192.168.0.150:7878/nph-mjpeg.cgi','','', 'Timeout', 20);
%set(cam1, 'ReturnedColorSpace', 'RGB');
imgFromCam1 = imread('C:\FFOutput\imgFromCamL_1.png');
 
% This is where and what my image will be saved
counter  = 1;
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\L_2\';
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
%cam2 = videoinput('winvideo', 2, 'RGB24_960x720');
%Cam_R = ipcam('http://192.168.0.198:80/mjpg/video.mjpg','root','0000', 'Timeout', 20);
%set(cam2, 'ReturnedColorSpace', 'RGB');
imgFromCam2 = imread('C:\FFOutput\imgFromCamR_1.png');
 
% This is where and what my image will be saved
counter  = 1;
baseDir  = 'C:\opencv\AR\Stereo camera\Fundamental matrix\R_2\';
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
% --- 校正
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
% --- 定位
function pushbutton6_Callback(hObject, eventdata, handles)
conn=database('ar_viewer','admin','winnie448t','com.mysql.jdbc.Driver','jdbc:mysql://140.127.34.215:3306/ar_viewer');
global LOC_OBJ_ NUM_OF_OBJ DIS_OF_OBJ_;
global img_L P1;
clc;

img_L = imread('C:\FFOutput\imgFromCamL_1.png');
img_R = imread('C:\FFOutput\imgFromCamR_1.png');

load('stereoParams.mat');

[imgRect1, imgRect2] = rectifyStereoImages(img_L,img_R,stereoParams);

imgGray1 = rgb2gray(imgRect1);
imgGray2 = rgb2gray(imgRect2);

disparityMap = disparity(imgGray1,imgGray2);

colormap jet
colorbar

point3D = reconstructScene(disparityMap, stereoParams);

z = point3D(:,:,3);
z(z < 0 | z > 8000) = NaN;
x = point3D(:,:,1);
x(x < -3000 | x > 3000) = NaN;
point3D(:,:,3) = z;
point3D(:,:,1) = x;

NUM_OF_OBJ=1;
imshow(img_L,'Parent',handles.axes2);
while NUM_OF_OBJ>=0
    objectRegion = round(getPosition(imrect)); % mouse click and drag
    confirm = questdlg('Is the area you choose correct?','Confirm','Yes','No','Yes');
    switch confirm
        case 'Yes'
            if NUM_OF_OBJ == 1
                object = insertShape(img_L, 'Rectangle', objectRegion, 'LineWidth',2,'Color','yellow');
                object2=insertText(object,objectRegion(1:2),NUM_OF_OBJ,'FontSize',12,'BoxColor','red','TextColor','white');
                imshow(object2,'Parent',handles.axes2); 

                centroids = [round(objectRegion(:, 1) + objectRegion(:,3 ) / 2),... 
                round(objectRegion(:, 2) + objectRegion(:, 4) / 2)];

                centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
                X = point3D(:, :, 1);
                Y = point3D(:, :, 2);
                Z = point3D(:, :, 3);
                centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];

                centroids3D = reshape(centroids3D,[1,3]); % reshape into 1 row, 3 column

                LOC_OBJ_{NUM_OF_OBJ}=num2str(centroids3D);
                str = deblank(num2str(centroids3D));
                S = regexp(str, '\s+', 'split');
                X = char(S{1});
                Y = char(S{2});
                Z = char(S{3});
                
                colnames = {'X','Y','Z','ID'};
                data={X,Y,Z,NUM_OF_OBJ};
                tablename = 'Location';
                insert(conn,tablename,colnames,data);

                dists = sqrt(sum(centroids3D .^ 2)) / 1000;

                labels = cell(1, numel(dists));
                for i = 1:numel(dists)
                    labels{i} = sprintf('Distance: %0.2f meters', dists(i));
                end
                DIS_OF_OBJ_{NUM_OF_OBJ}=labels;
                colnames_2 = {'distance','ID'};
                data_2={dists,NUM_OF_OBJ};
                tablename_2 = 'Distance';
                insert(conn,tablename_2,colnames_2,data_2);
                
                pts1 = detectMinEigenFeatures(imgGray1, 'ROI', objectRegion,'MinQuality',0.1);
                % Create the points tracker.
                tracker = vision.PointTracker('MaxBidirectionalError',1,'NumPyramidLevels',5);

                % Initialize the point tracker.
                pts1 = pts1.Location;
                initialize(tracker,pts1,imgGray1);

                % Track the points.
                [pts2, validIdx] = step(tracker, imgGray2);

                % Compute 3-D locations of matching points in two images
                worldPoints = triangulate(pts1,pts2,stereoParams);

                % Compute R and t
                [R1,t1] = extrinsics(pts1,worldPoints,stereoParams.CameraParameters1);
                [R2,t2] = extrinsics(pts2,worldPoints,stereoParams.CameraParameters2);

                % Compute camera projection matrix
                P1 = cameraMatrix(stereoParams.CameraParameters1, R1, t1);
                P2 = cameraMatrix(stereoParams.CameraParameters2, R2, t2);
                makecube = questdlg('You have already make a cube.');
                
            else
                object2 = insertShape(object2, 'Rectangle', objectRegion, 'LineWidth',2,'Color','yellow');
                object2=insertText(object2,objectRegion(1:2),NUM_OF_OBJ,'FontSize',12,'BoxColor','red','TextColor','white');
                imshow(object2,'Parent',handles.axes2);
                centroids = [round(objectRegion(:, 1) + objectRegion(:,3 ) / 2),... 
                round(objectRegion(:, 2) + objectRegion(:, 4) / 2)];

                centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), centroids(:, 1));
                X = point3D(:, :, 1);
                Y = point3D(:, :, 2);
                Z = point3D(:, :, 3);
                centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];

                centroids3D = reshape(centroids3D,[1,3]);

                LOC_OBJ_{NUM_OF_OBJ}=num2str(centroids3D);
                str = deblank(num2str(centroids3D));
                S = regexp(str, '\s+', 'split');
                X = char(S{1});
                Y = char(S{2});
                Z = char(S{3});
                
                colnames = {'X','Y','Z','ID'};
                data={X,Y,Z,NUM_OF_OBJ};
                tablename = 'Location';
                insert(conn,tablename,colnames,data);
                
                dists = sqrt(sum(centroids3D .^ 2)) / 1000;
                
                labels = cell(1, numel(dists));
                for i = 1:numel(dists)
                    labels{i} = sprintf('Distance: %0.2f meters', dists(i));
                end
                DIS_OF_OBJ_{NUM_OF_OBJ}=labels;
                colnames_2 = {'distance','ID'};
                data_2={dists,NUM_OF_OBJ};
                tablename_2 = 'Distance';
                insert(conn,tablename_2,colnames_2,data_2);
                
                pts1 = detectMinEigenFeatures(imgGray1, 'ROI', objectRegion,'MinQuality',0.1);
                % Create the points tracker.
                tracker = vision.PointTracker('MaxBidirectionalError',1,'NumPyramidLevels',5);

                % Initialize the point tracker.
                pts1 = pts1.Location;
                initialize(tracker,pts1,imgGray1);

                % Track the points.
                [pts2, validIdx] = step(tracker, imgGray2);

                % Compute 3-D locations of matching points in two images
                worldPoints = triangulate(pts1,pts2,stereoParams);

                % Compute R and t
                [R1,t1] = extrinsics(pts1,worldPoints,stereoParams.CameraParameters1);
                [R2,t2] = extrinsics(pts2,worldPoints,stereoParams.CameraParameters2);

                % Compute camera projection matrix
                P1 = cameraMatrix(stereoParams.CameraParameters1, R1, t1);
                P2 = cameraMatrix(stereoParams.CameraParameters2, R2, t2);
                makecube = questdlg('You have already make a cube.');
            end
            NUM_OF_OBJ=NUM_OF_OBJ+1;
 
        case 'No'
            if NUM_OF_OBJ == 1
                NUM_OF_OBJ=1;
            end
    end
end

% --- 導覽
function pushbutton8_Callback(hObject, eventdata, handles)
set(handles.text10,'Visible','On');
set(handles.text11,'Visible','On');
set(handles.text12,'Visible','On');
set(handles.text9,'Visible','On');
set(handles.text8,'Visible','On');
set(handles.popupmenu2,'Visible','On');
set(handles.pushbutton13,'Visible','On');
conn=database('ar_viewer','admin','winnie448t','com.mysql.jdbc.Driver','jdbc:mysql://140.127.34.215:3306/ar_viewer');
sqlquery = 'SELECT MAX(ID) FROM location';
curs = exec(conn,sqlquery);
curs = fetch(curs);
max = cell2mat(curs.Data);

for i=1:max
    list{i}=i;
end
set(handles.popupmenu2,'string',list);



    


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- 展示資料
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_L P1 ;
contents = get(handles.popupmenu2,'String'); 
popupmenu4value = contents{get(handles.popupmenu2,'Value')};
conn=database('ar_viewer','admin','winnie448t','com.mysql.jdbc.Driver','jdbc:mysql://140.127.34.215:3306/ar_viewer');

switch popupmenu4value
    case '1'
                sqlquery = 'SELECT X,Y,Z FROM location WHERE ID=1';
                sqlquery_dis = 'SELECT distance FROM distance WHERE ID=1';
                curs = exec(conn,sqlquery);
                curs_dis = exec(conn,sqlquery_dis);
                curs = fetch(curs);
                curs_dis = fetch(curs_dis);
                data = curs.Data;
                data_dis = curs_dis.Data;
                set(handles.text9,'String',data);
                set(handles.text10,'String',data_dis);
                a = cellfun(@str2num,data);
                 %==========================================================================
                % Find the 8 vertices of cube around 3-D point
                %==========================================================================
                w = 200; % Full-width of each side of cube.
                h = w/2; % For readability.

                figure('Visible','Off');
                % Face #1 of 6.
                cube.f1 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[-h -h -h -h], ...
                    'FaceColor', 'w', 'FaceAlpha', 0.1 ); % bottom
                daspect( [1 1 1] )  % 1:1:1 aspect ratio.
                hold on

                % Face #2 of 6.
                cube.f2 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[ h  h  h  h], ...
                    'FaceColor', 'r', 'FaceAlpha', 0.7 ); % top

                % Face #3 of 6.
                cube.f3 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[ h  h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'y', 'FaceAlpha', 0.7 );

                % Face #4 of 6.
                cube.f4 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h -h -h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'b', 'FaceAlpha', 0.7 );

                % Face #5 of 6.
                cube.f5 = patch( 'XData', a(1)+[ 2*h  2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'k', 'FaceAlpha', 0.7 );

                % Face #6 of 6.
                cube.f6 = patch( 'XData', a(1)+[-2*h -2*h -2*h -2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'g', 'FaceAlpha', 0.7 );

                % Display 3-D point in the middle of cube.
                scatter3(a(1), a(2), a(3),'or', 'filled', 'SizeData', 50 )
                xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
                view(3)
                title('A cube around 3-D point')
                grid on

                % A cube. A cube is defined by 8 vertices that form 6 faces. Eeach of the 6
                % faces has 4 vertices.
                face1to4 = get(cube.f1,'Vertices'); % 4 vertices from bottom
                face5to8 = get(cube.f2,'Vertices'); % 4 vertices from top
                vertices = cat(1,face1to4,face5to8); % 8 vertices from bottom and top

                obj.vertices = [vertices(1,:); % first vertice
                     vertices(2,:);
                     vertices(3,:);
                     vertices(4,:);
                     vertices(5,:);
                     vertices(6,:);
                     vertices(7,:);
                     vertices(8,:);];
                 faces = ...
                    [1 4 3 2; % first face
                    5 8 7 6;
                    5 1 2 6;
                    8 4 3 7;
                    6 2 3 7;
                    5 1 4 8;];
                %==========================================================================
                % 3D-to-2D projection (Cam 1)
                %==========================================================================
                % Homogeneous coordinates
                vert = [obj.vertices,ones(size(obj.vertices,1),1)];
               
                % Transform vertices from 3D world into 2D image plane
                vert12D = vert*P1;

                % Homogenous coordinates to (x,y)      
                for i1 = 1:size(vert12D,1)
                   vertImage1(i1,1:2) =  vert12D(i1,1:2)/vert12D(i1,3);
                end 
                % Draw patch/display a cube                                    
                figure();
                imshow(img_L);
                hold on
                obj.patch1 = patch('Faces',faces,'Vertices',vertImage1,'EdgeColor','b','facecolor','none','LineWidth',1);
                title('AR viewing 3D-to-2D projection (Cam 1)');

                % Convert 3D point to 2D, then display in the middle of cube
                p1 = [a,ones(size(a,1),1)];
                p12D = p1*P1;
                for j1 = 1:size(p12D,1)
                    vImage1(j1,1:2) = p12D(j1,1:2)/p12D(j1,3);
                end
                    scatter(vImage1(1),vImage1(2),'MarkerFaceColor','r')
    case '2'
                sqlquery = 'SELECT X,Y,Z FROM location WHERE ID=2';
                sqlquery_dis = 'SELECT distance FROM distance WHERE ID=2';
                curs = exec(conn,sqlquery);
                curs_dis = exec(conn,sqlquery_dis);
                curs = fetch(curs);
                curs_dis = fetch(curs_dis);
                data = curs.Data;
                data_dis = curs_dis.Data;
                set(handles.text9,'String',data);
                set(handles.text10,'String',data_dis);
                a = cellfun(@str2num,data);
                %==========================================================================
                % Find the 8 vertices of cube around 3-D point
                %==========================================================================
                w = 200; % Full-width of each side of cube.
                h = w/2; % For readability.

                figure('Visible','Off');
                % Face #1 of 6.
                cube.f1 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[-h -h -h -h], ...
                    'FaceColor', 'w', 'FaceAlpha', 0.1 ); % bottom
                daspect( [1 1 1] )  % 1:1:1 aspect ratio.
                hold on

                % Face #2 of 6.
                cube.f2 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[ h  h  h  h], ...
                    'FaceColor', 'r', 'FaceAlpha', 0.7 ); % top

                % Face #3 of 6.
                cube.f3 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[ h  h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'y', 'FaceAlpha', 0.7 );

                % Face #4 of 6.
                cube.f4 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h -h -h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'b', 'FaceAlpha', 0.7 );

                % Face #5 of 6.
                cube.f5 = patch( 'XData', a(1)+[ 2*h  2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'k', 'FaceAlpha', 0.7 );

                % Face #6 of 6.
                cube.f6 = patch( 'XData', a(1)+[-2*h -2*h -2*h -2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'g', 'FaceAlpha', 0.7 );

                % Display 3-D point in the middle of cube.
                scatter3(a(1), a(2), a(3),'or', 'filled', 'SizeData', 50 )
                xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
                view(3)
                title('A cube around 3-D point')
                grid on

                % A cube. A cube is defined by 8 vertices that form 6 faces. Eeach of the 6
                % faces has 4 vertices.
                face1to4 = get(cube.f1,'Vertices'); % 4 vertices from bottom
                face5to8 = get(cube.f2,'Vertices'); % 4 vertices from top
                vertices = cat(1,face1to4,face5to8); % 8 vertices from bottom and top

                obj.vertices = [vertices(1,:); % first vertice
                     vertices(2,:);
                     vertices(3,:);
                     vertices(4,:);
                     vertices(5,:);
                     vertices(6,:);
                     vertices(7,:);
                     vertices(8,:);];
                 faces = ...
                    [1 4 3 2; % first face
                    5 8 7 6;
                    5 1 2 6;
                    8 4 3 7;
                    6 2 3 7;
                    5 1 4 8;];
                %==========================================================================
                % 3D-to-2D projection (Cam 1)
                %==========================================================================
                % Homogeneous coordinates
                vert = [obj.vertices,ones(size(obj.vertices,1),1)];
               
                % Transform vertices from 3D world into 2D image plane
                vert12D = vert*P1;

                % Homogenous coordinates to (x,y)      
                for i1 = 1:size(vert12D,1)
                   vertImage1(i1,1:2) =  vert12D(i1,1:2)/vert12D(i1,3);
                end 
                % Draw patch/display a cube                                    
                figure();
                imshow(img_L);
                hold on
                obj.patch1 = patch('Faces',faces,'Vertices',vertImage1,'EdgeColor','b','facecolor','none','LineWidth',1);
                title('AR viewing 3D-to-2D projection (Cam 1)');

                % Convert 3D point to 2D, then display in the middle of cube
                a = cellfun(@str2num,data);
                p1 = [a,ones(size(a,1),1)];
                p12D = p1*P1;
                for j1 = 1:size(p12D,1)
                    vImage1(j1,1:2) = p12D(j1,1:2)/p12D(j1,3);
                end
                    scatter(vImage1(1),vImage1(2),'MarkerFaceColor','r')
    case '3'
                sqlquery = 'SELECT X,Y,Z FROM location WHERE ID=3';
                sqlquery_dis = 'SELECT distance FROM distance WHERE ID=3';
                curs = exec(conn,sqlquery);
                curs_dis = exec(conn,sqlquery_dis);
                curs = fetch(curs);
                curs_dis = fetch(curs_dis);
                data = curs.Data;
                data_dis = curs_dis.Data;
                set(handles.text9,'String',data);
                set(handles.text10,'String',data_dis);
                a = cellfun(@str2num,data);
                %==========================================================================
                % Find the 8 vertices of cube around 3-D point
                %==========================================================================
                w = 200; % Full-width of each side of cube.
                h = w/2; % For readability.

                figure('Visible','Off');
                % Face #1 of 6.
                cube.f1 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[-h -h -h -h], ...
                    'FaceColor', 'w', 'FaceAlpha', 0.1 ); % bottom
                daspect( [1 1 1] )  % 1:1:1 aspect ratio.
                hold on

                % Face #2 of 6.
                cube.f2 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[ h  h  h  h], ...
                    'FaceColor', 'r', 'FaceAlpha', 0.7 ); % top

                % Face #3 of 6.
                cube.f3 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[ h  h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'y', 'FaceAlpha', 0.7 );

                % Face #4 of 6.
                cube.f4 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h -h -h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'b', 'FaceAlpha', 0.7 );

                % Face #5 of 6.
                cube.f5 = patch( 'XData', a(1)+[ 2*h  2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'k', 'FaceAlpha', 0.7 );

                % Face #6 of 6.
                cube.f6 = patch( 'XData', a(1)+[-2*h -2*h -2*h -2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'g', 'FaceAlpha', 0.7 );

                % Display 3-D point in the middle of cube.
                scatter3(a(1), a(2), a(3),'or', 'filled', 'SizeData', 50 )
                xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
                view(3)
                title('A cube around 3-D point')
                grid on

                % A cube. A cube is defined by 8 vertices that form 6 faces. Eeach of the 6
                % faces has 4 vertices.
                face1to4 = get(cube.f1,'Vertices'); % 4 vertices from bottom
                face5to8 = get(cube.f2,'Vertices'); % 4 vertices from top
                vertices = cat(1,face1to4,face5to8); % 8 vertices from bottom and top

                obj.vertices = [vertices(1,:); % first vertice
                     vertices(2,:);
                     vertices(3,:);
                     vertices(4,:);
                     vertices(5,:);
                     vertices(6,:);
                     vertices(7,:);
                     vertices(8,:);];
                 faces = ...
                    [1 4 3 2; % first face
                    5 8 7 6;
                    5 1 2 6;
                    8 4 3 7;
                    6 2 3 7;
                    5 1 4 8;];
                %==========================================================================
                % 3D-to-2D projection (Cam 1)
                %==========================================================================
                % Homogeneous coordinates
                vert = [obj.vertices,ones(size(obj.vertices,1),1)];
               
                % Transform vertices from 3D world into 2D image plane
                vert12D = vert*P1;

                % Homogenous coordinates to (x,y)      
                for i1 = 1:size(vert12D,1)
                   vertImage1(i1,1:2) =  vert12D(i1,1:2)/vert12D(i1,3);
                end 
                % Draw patch/display a cube                                    
                figure();
                imshow(img_L);
                hold on
                obj.patch1 = patch('Faces',faces,'Vertices',vertImage1,'EdgeColor','b','facecolor','none','LineWidth',1);
                title('AR viewing 3D-to-2D projection (Cam 1)');

                % Convert 3D point to 2D, then display in the middle of cube
                a = cellfun(@str2num,data);
                p1 = [a,ones(size(a,1),1)];
                p12D = p1*P1;
                for j1 = 1:size(p12D,1)
                    vImage1(j1,1:2) = p12D(j1,1:2)/p12D(j1,3);
                end
                    scatter(vImage1(1),vImage1(2),'MarkerFaceColor','r')
                    
                    
                    
                     case '1'
                sqlquery = 'SELECT X,Y,Z FROM location WHERE ID=1';
                sqlquery_dis = 'SELECT distance FROM distance WHERE ID=1';
                curs = exec(conn,sqlquery);
                curs_dis = exec(conn,sqlquery_dis);
                curs = fetch(curs);
                curs_dis = fetch(curs_dis);
                data = curs.Data;
                data_dis = curs_dis.Data;
                set(handles.text9,'String',data);
                set(handles.text10,'String',data_dis);
                a = cellfun(@str2num,data);
                 %==========================================================================
                % Find the 8 vertices of cube around 3-D point
                %==========================================================================
                w = 200; % Full-width of each side of cube.
                h = w/2; % For readability.

                figure('Visible','Off');
                % Face #1 of 6.
                cube.f1 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[-h -h -h -h], ...
                    'FaceColor', 'w', 'FaceAlpha', 0.1 ); % bottom
                daspect( [1 1 1] )  % 1:1:1 aspect ratio.
                hold on

                % Face #2 of 6.
                cube.f2 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h  h  h -h], ...
                    'ZData', a(3)+2*[ h  h  h  h], ...
                    'FaceColor', 'r', 'FaceAlpha', 0.7 ); % top

                % Face #3 of 6.
                cube.f3 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[ h  h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'y', 'FaceAlpha', 0.7 );

                % Face #4 of 6.
                cube.f4 = patch( 'XData', a(1)+[-2*h -2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h -h -h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'b', 'FaceAlpha', 0.7 );

                % Face #5 of 6.
                cube.f5 = patch( 'XData', a(1)+[ 2*h  2*h  2*h  2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'k', 'FaceAlpha', 0.7 );

                % Face #6 of 6.
                cube.f6 = patch( 'XData', a(1)+[-2*h -2*h -2*h -2*h], ...
                    'YData', a(2)+2*[-h -h  h  h], ...
                    'ZData', a(3)+2*[-h  h  h -h], ...
                    'FaceColor', 'g', 'FaceAlpha', 0.7 );

                % Display 3-D point in the middle of cube.
                scatter3(a(1), a(2), a(3),'or', 'filled', 'SizeData', 50 )
                xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
                view(3)
                title('A cube around 3-D point')
                grid on

                % A cube. A cube is defined by 8 vertices that form 6 faces. Eeach of the 6
                % faces has 4 vertices.
                face1to4 = get(cube.f1,'Vertices'); % 4 vertices from bottom
                face5to8 = get(cube.f2,'Vertices'); % 4 vertices from top
                vertices = cat(1,face1to4,face5to8); % 8 vertices from bottom and top

                obj.vertices = [vertices(1,:); % first vertice
                     vertices(2,:);
                     vertices(3,:);
                     vertices(4,:);
                     vertices(5,:);
                     vertices(6,:);
                     vertices(7,:);
                     vertices(8,:);];
                 faces = ...
                    [1 4 3 2; % first face
                    5 8 7 6;
                    5 1 2 6;
                    8 4 3 7;
                    6 2 3 7;
                    5 1 4 8;];
                %==========================================================================
                % 3D-to-2D projection (Cam 1)
                %==========================================================================
                % Homogeneous coordinates
                vert = [obj.vertices,ones(size(obj.vertices,1),1)];
               
                % Transform vertices from 3D world into 2D image plane
                vert12D = vert*P1;

                % Homogenous coordinates to (x,y)      
                for i1 = 1:size(vert12D,1)
                   vertImage1(i1,1:2) =  vert12D(i1,1:2)/vert12D(i1,3);
                end 
                % Draw patch/display a cube                                    
                figure();
                imshow(img_L);
                hold on
                obj.patch1 = patch('Faces',faces,'Vertices',vertImage1,'EdgeColor','b','facecolor','none','LineWidth',1);
                title('AR viewing 3D-to-2D projection (Cam 1)');

                % Convert 3D point to 2D, then display in the middle of cube
                p1 = [a,ones(size(a,1),1)];
                p12D = p1*P1;
                for j1 = 1:size(p12D,1)
                    vImage1(j1,1:2) = p12D(j1,1:2)/p12D(j1,3);
                end
                    scatter(vImage1(1),vImage1(2),'MarkerFaceColor','r')
end
