%==========================================================================
% This test shows how to defind object of interest with a calibrated stereo
% camera and determine distance from the cameras.
%==========================================================================
clear;
close all;
clc;
%==========================================================================
% Create the webcam object.
%==========================================================================
cam1 = videoinput('winvideo', 1, 'RGB24_640x480'); % 1 = 1st camera
cam2 = videoinput('winvideo', 2, 'RGB24_640x480'); % 2 = 2nd camera

% Capture one frame to get its size.
img1 = getsnapshot(cam1); frameSize1 = size(img1);
img2 = getsnapshot(cam2); frameSize2 = size(img2);

%==========================================================================
% Load precomputed stereoParameters object. 
%==========================================================================
% Load the |stereoParameters| object, which is the result of calibrating
% the camera.
load('stereoParams.mat');

%==========================================================================
% Rectify images
%==========================================================================
% The images from the first and the second cameras must be rectified in order
% to compute disparity and reconstruct the 3-D scene. Rectified images have
% horizontal epipolar lines, and are row-aligned. This simplifies the
% computation of disparity by reducing the search space for matching points
% to one dimension. Rectified images can also be combined into an anaglyph,
% which can be viewed using the stereo red-cyan glasses to see the 3-D
% effect.
[imgRect1, imgRect2] = rectifyStereoImages(img1,img2,stereoParams,'OutputView','full');

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

% Limit the range of Z and X for display.
z = point3D(:,:,3);
z(z < 0 | z > 8000) = NaN;
x = point3D(:,:,1);
x(x < -3000 | x > 3000) = NaN;
point3D(:,:,3) = z;
point3D(:,:,1) = x;

figure;
pcshow(point3D, imgRect1, 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
set(gca, 'CameraViewAngle',10, 'CameraUpVector',[0 -1 0],...
    'CameraPosition',[16500 -13852 -49597], 'DataAspectRatio',[1 1 1]);
title('Reconstructed 3-D Scene');

%==========================================================================
% Detect object of interest
%==========================================================================
% Define object region using a mouse in the first image
figure; imshow(imgGray1);
objectRegion = round(getPosition(imrect)); % mouse click and drag

% Show object region with a bounding box in the first camera.
object = insertShape(imgGray1, 'Rectangle', objectRegion, 'LineWidth',4,'Color','red');
figure; imshow(object); title('Red Box Shows Object Region');
%{ 
% Detect interest points in the object region.
pts1 = detectMinEigenFeatures(imgGray1, 'ROI', objectRegion,'MinQuality',0.1);

% Visualize the detected points.
figure; imshow(imgGray1,'InitialMagnification',50);
title('150 Strongest corners from the first cam');
hold on
plot(selectStrongest(pts1,150));

% Create the points tracker.
tracker = vision.PointTracker('MaxBidirectionalError',1,'NumPyramidLevels',5);

% Initialize the point tracker.
pts1 = pts1.Location;
initialize(tracker,pts1,imgGray1);

% Track the points.
[pts2, validIdx] = step(tracker, imgGray2);
matchedPts1 = pts1(validIdx, :);
matchedPts2 = pts2(validIdx, :);

% Visualize correspondences
figure; showMatchedFeatures(imgGray1, imgGray2, matchedPts1, matchedPts2, 'montage');
title('Tracked Features');
%}
%==========================================================================
% Find the 3-D world coordinates of the centroid of object and compute the
% distance from the centroid to the camera in meters.
%==========================================================================

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
disp3D = insertObjectAnnotation(imgRect1,'rectangle', objectRegion,...
    ['XYZ: ',num2str(reshapeCentroids3D)], 'FontSize', 18);
figure;
imshow(disp3D);title('Detected object with 3-D world coordinates (scale in millimeter)');

% Find the distances from the camera in meters.
dists = sqrt(sum(centroids3D .^ 2)) / 1000;

% Display the detected object and distances.
labels = cell(1, numel(dists));
for i = 1:numel(dists)
    labels{i} = sprintf('Distance: %0.2f meters', dists(i));
end
figure;
imshow(insertObjectAnnotation(imgRect1, 'rectangle', objectRegion, labels,'FontSize', 18));
title('Detected Object with distance (scale in meter)');

rgb = insertMarker(imgRect1,centroids,'+','size',10);
figure; imshow(rgb);


