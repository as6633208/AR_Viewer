clear;      % clears all variables in the MatLab workspace
clc;        % clears the command window
close all;  % closes all figure windows

%==========================================================================
% Create the webcam object.
%==========================================================================
%cam1 = videoinput('winvideo', 1, 'RGB24_640x480'); 
%cam2 = videoinput('winvideo', 2, 'RGB24_640x480');
img_L = imread('C:\FFOutput\imgFromCamL_1.png');
img_R = imread('C:\FFOutput\imgFromCamR_1.png');

% Capture one frame to get its size.
%img1 = getsnapshot(cam1);frameSize1 = size(img1); 
%img2 = getsnapshot(cam2);frameSize2 = size(img2); 

%==========================================================================
% Load precomputed stereoParameters object. 
%==========================================================================
% Load the |stereoParameters| object, which is the result of calibrating
% the camera.
load('stereoParams.mat');

%==========================================================================
% Rectify images
%==========================================================================
% The images from the first and the second cameras must be rectified in
% order to compute disparity and reconstruct the 3-D scene. Rectified
% images have horizontal epipolar lines, and are row-alligned. This
% simplifies the computation of disparity by reducing the search space for
% matching points to one dimension. Rectified images can also be comined
% into an anaglyph, which can be viewd using the stereo red-cyan glasses to
% see the 3-D effect.
[imgRect1, imgRect2] = rectifyStereoImages(img_L,img_R,stereoParams);

figure;
imshow(stereoAnaglyph(imgRect1, imgRect2)); title('Rectified Image');

%==========================================================================
% Compute disparity
%==========================================================================
% In rectified stereo images any pair of corresponding points are located
% on the same pixel row.  For each pixel in the left image compute the
% distance to the correspoinding pixel in the right image. This distance is
% called the disparity, and it is proportional to the distance of the
% corresponding world point from the camera.
imgGray1 = rgb2gray(imgRect1);
imgGray2 = rgb2gray(imgRect2);

disparityMap = disparity(imgGray1,imgGray2);
figure;
imshow(disparityMap, [0, 64]); title('Disparity Map');
colormap jet
colorbar

%==========================================================================
% Reconstruct the 3-D Scene
%==========================================================================
% Reconstruct the 3-D world coordinates of points corresponding to each
% pixel from the disparity map.
point3D = reconstructScene(disparityMap, stereoParams);

figure;
pcshow(point3D, imgRect1, 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
set(gca, 'CameraViewAngle',10, 'CameraUpVector',[0 -1 0],...
    'CameraPosition',[16500 -13852 -49597], 'DataAspectRatio',[1 1 1]);
title('Reconstructed 3-D Scene');

%==========================================================================
% Select the object region using a mouse 
%==========================================================================
% Select the object region using a mouse in the first image
figure; imshow(imgGray1);
objectRegion = round(getPosition(imrect)); % mouse click and drag
title('Select the Object Region Using a Mouse');

% Show object region with a bounding box in the first camera.
object1 = insertShape(imgGray1, 'Rectangle', objectRegion, 'LineWidth',...
    4,'Color','red');
figure; imshow(object1); title('Red Box Shows Object Region');
hold on

%==========================================================================
% Find the 3-D world coordinates of the centroid of object
%==========================================================================
% Find the centroids of object
centroids = [round(objectRegion(:, 1) + objectRegion(:,3 ) / 2),... 
    round(objectRegion(:, 2) + objectRegion(:, 4) / 2)];

% Display the centroids of object
rgb = insertMarker(imgRect1,centroids,'+','size',10); 
figure; imshow(rgb);
title('The centroids of object');

% Find the 3-D world coordinates of the centroids
centroidsIdx = sub2ind(size(disparityMap), centroids(:, 2), ...
    centroids(:, 1));
X = point3D(:, :, 1);
Y = point3D(:, :, 2);
Z = point3D(:, :, 3);
centroids3D = [X(centroidsIdx)'; Y(centroidsIdx)'; Z(centroidsIdx)'];
centroids3D = reshape(centroids3D,[1,3]); % reshape into 1 row, 3 column

% Display detected object with 3-D world coordinates
disp3D = insertObjectAnnotation(imgRect1,'rectangle', objectRegion,...
    ['XYZ: ',num2str(centroids3D)], 'FontSize', 18);
figure;
imshow(disp3D);
title('Detected object with 3-D world coordinates (scale in millimeter)');

%==========================================================================
% Compute the distance from the centroid to the camera in meters
%==========================================================================
% Find the distances from the camera in meters.
dists = sqrt(sum(centroids3D .^ 2)) / 1000;

% Display the detected object and distances.
labels = cell(1, numel(dists));
for i1 = 1:numel(dists)
    labels{i1} = sprintf('Distance: %0.2f meters', dists(i1));
end
figure;
imshow(insertObjectAnnotation(imgRect1, 'rectangle', objectRegion, ...
    labels,'FontSize', 18));
title('Detected object with distance (scale in meter)');

%==========================================================================
% Compute camera projection matrix
%==========================================================================
%{
Syntax
camMat = cameraMatrix(cameraParams,rotationMatrix,translationVector)

Input:
    1. cameraParams - camera parameters

    2. rotationMatrix - rotation of camera
    Rotation of camera, specified as a 3-by-3 matrix. 

    3. translationVector - translation of camera
    Translation of camera, specified as a 1-by-3 element vector. The 
    translation vector describes the transformation from the world 
    coordinates to the camera coordinates. 

Output:
    1. camProjMat - camera projection matrix
    Camera projection matrix, returned as a 4-by-3 matrix. The matrix
    contains the 3-D world points in homogenous coordinates that are
    projected into the image. The function computes camProjMat as follows:
        
        camProjMat = [rotationMatrix; translationVector] ?K.
        K: the intrinsic matrix

Then, using the camera matrix and homogeneous coordinates, you can project 
a world point onto the image.

        w ?[x,y,1] = [X,Y,Z,1] ?camMatrix.
(X,Y,Z): world coordinates of a point
(x,y): coordinates of the corresponding image point
w: arbitrary scale factor
%}
% Detect interest points in the object region.
pts1 = detectMinEigenFeatures(imgGray1, 'ROI', objectRegion,...
    'MinQuality',0.1);

% Visualize the detected points.
figure; imshow(imgGray1,'InitialMagnification',50);
title('150 Strongest corners from the first cam');
hold on
plot(selectStrongest(pts1,150));

% Create the points tracker.
tracker = vision.PointTracker('MaxBidirectionalError',1,...
    'NumPyramidLevels',5);

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
%==========================================================================
% Find the 8 vertices of cube around 3-D point
%==========================================================================
w = 250; % Full-width of each face of cube.
h = w/2; % For readability.

figure;
% Face #1 of 6. Bottom
cube.f1 = patch( 'XData', centroids3D(1)+[-h -h  h  h], ...
    'YData', centroids3D(2)+[-h  h  h -h], ...
    'ZData', centroids3D(3)+[-h -h -h -h], ...
    'FaceColor', 'w', 'FaceAlpha', 0.1 ); 
daspect( [1 1 1] )  % 1:1:1 aspect ratio.
hold on

% Face #2 of 6. Top
cube.f2 = patch( 'XData', centroids3D(1)+[-h -h  h  h], ...
    'YData', centroids3D(2)+[-h  h  h -h], ...
    'ZData', centroids3D(3)+[ h  h  h  h], ...
    'FaceColor', 'r', 'FaceAlpha', 0.7 ); 

% Face #3 of 6.
cube.f3 = patch( 'XData', centroids3D(1)+[-h -h  h  h], ...
    'YData', centroids3D(2)+[ h  h  h  h], ...
    'ZData', centroids3D(3)+[-h  h  h -h], ...
    'FaceColor', 'y', 'FaceAlpha', 0.7 );

% Face #4 of 6.
cube.f4 = patch( 'XData', centroids3D(1)+[-h -h  h  h], ...
    'YData', centroids3D(2)+[-h -h -h -h], ...
    'ZData', centroids3D(3)+[-h  h  h -h], ...
    'FaceColor', 'b', 'FaceAlpha', 0.7 );

% Face #5 of 6.
cube.f5 = patch( 'XData', centroids3D(1)+[ h  h  h  h], ...
    'YData', centroids3D(2)+[-h -h  h  h], ...
    'ZData', centroids3D(3)+[-h  h  h -h], ...
    'FaceColor', 'k', 'FaceAlpha', 0.7 );

% Face #6 of 6.
cube.f6 = patch( 'XData', centroids3D(1)+[-h -h -h -h], ...
    'YData', centroids3D(2)+[-h -h  h  h], ...
    'ZData', centroids3D(3)+[-h  h  h -h], ...
    'FaceColor', 'g', 'FaceAlpha', 0.7 );

% Display 3-D point in the middle of cube.
scatter3(centroids3D(1), centroids3D(2), centroids3D(3), ...
    'or', 'filled', 'SizeData', 50 )
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

obj.faces = ...
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
figure;
imshow(img_L);
hold on
obj.patch1 = patch('Faces',obj.faces,'Vertices',vertImage1,...
    'EdgeColor','r','facecolor','none','LineWidth',1);
title('AR viewing 3D-to-2D projection (Cam 1)');

% Convert 3D point to 2D, then display in the middle of cube
p1 = [centroids3D,ones(size(centroids3D,1),1)];
p12D = p1*P1;
for j1 = 1:size(p12D,1)
    vImage1(j1,1:2) = p12D(j1,1:2)/p12D(j1,3);
end
scatter(vImage1(1),vImage1(2),'MarkerFaceColor','r')

% Display 8 vertices
%x1 = vertImage1(:,1); % extract the frist column
%y1 = vertImage1(:,end); % extract the last/second column
%scatter(x1,y1,'MarkerFaceColor','b')

%{
% Display the coordinate of the 8 vertices
x = double(x); % convert to double precision
y = double(y);
for k =1:numel(x)
    text(x(k),y(k),['(' num2str(x(k)) ',' num2str(y(k)) ')'],'Color','red')
end
%}

%==========================================================================
% 3D-to-2D projection (Cam 2)
%==========================================================================
% Homogeneous coordinates
%vert = [obj.vertices,ones(size(obj.vertices,1),1)];

% Transform vertices from 3D world into 2D image plane
vert22D = vert*P2;

% Homogenous coordinates to (x,y)
for i2 = 1:size(vert22D,1)
    vertImage2(i2,1:2) =  vert22D(i2,1:2)/vert22D(i2,3);
end 

% Draw patch/display a cube
figure;
imshow(img_R);
hold on
obj.patch2 = patch('Faces',obj.faces,'Vertices',vertImage2,...
    'EdgeColor','r','facecolor','none','LineWidth',1);
title('AR viewing 3D-to-2D projection (Cam 2)');

% Convert 3D point to 2D, then display in the middle of cube
p2 = [centroids3D,ones(size(centroids3D,1),1)];
p22D = p2*P2;
for j2 = 1:size(p22D,1)
    vImage2(j2,1:2) = p22D(j2,1:2)/p22D(j2,3);
end
scatter(vImage2(1),vImage2(2),'MarkerFaceColor','r')

% Display 8 vertices
%x2 = vertImage2(:,1); % extract the frist column
%y2 = vertImage2(:,end); % extract the last/second column
%scatter(x2,y2,'MarkerFaceColor','b')

%{
% Display the coordinate of the 8 vertices
x = double(x); % convert to double precision
y = double(y);
for k =1:numel(x)
    text(x(k),y(k),['(' num2str(x(k)) ',' num2str(y(k)) ')'],'Color','red')
end
%}
%==========================================================================
% Click at the object/cube, show some data on the object
%==========================================================================

c = plot(vImage2(1),vImage2(2),...
    'bo',...
    'Color', 'none',...
    'MarkerFaceColor', 'b',...
    'MarkerSize', 10); 
hold on

[xin,yin,rin]=MagnetGInput(c);

hold on;
plot(xin,yin,'ro', 'LineWidth',3,'MarkerSize',8);
text(xin+10,yin+.2,{'This is object A',...
    '.......................',...
    '.......................'},'color','r','FontSize',14);

plot(xin,yin,'ro', 'LineWidth',3,'MarkerSize',8);

[x,y] = ds2nfu(xin,yin); % convert to normalized fiure units
%{
annotation('textbox',...
    [x y 0.00 0.00],...
    'String',{'This is object A',...
    '.......................',...
    '.......................'},...
    'FitBoxToText','on',...
    'FontSize',14,...
    'EdgeColor','none',...
    'FontName','Arial',...
    'BackgroundColor','w',...
    'FaceAlpha',.7,...
    'Color','k');
    
%}




