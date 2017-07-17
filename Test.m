%==========================================================================
% Setup
%==========================================================================
% Create the webcam object.
cam1 = videoinput('winvideo', 1, 'RGB24_640x480');
cam2 = videoinput('winvideo', 2, 'RGB24_640x480'); 

% Capture one frame to get its size.
img1 = getsnapshot(cam1); frameSize1 = size(img1);
img2 = getsnapshot(cam2); frameSize2 = size(img2);

% Load precomputed stereo parameters.
load ('stereoParams.mat');

% Remove lens distortion.
img1 = undistortImage(img1,stereoParams.CameraParameters1);
img2 = undistortImage(img2,stereoParams.CameraParameters2);
figure(1); imshow(img1); hold on
figure(2); imshow(img2); hold on

% Convert from RGB to GRAY.
gr1 = rgb2gray(img1);
gr2 = rgb2gray(img2);

%==========================================================================
% Find point correspondences between camera1 and camera2
%==========================================================================
% Detect feature points.
points1 = detectMinEigenFeatures(gr1); 
points2 = detectMinEigenFeatures(gr2); 

% Extract the neighborhood features.
[features1, valid_points1] = extractFeatures(gr1, points1);
[features2, valid_points2] = extractFeatures(gr2, points2);

% Use the sum of absolute differences (SAD) metric to determine indices of
% matching features.
indexPairs = matchFeatures(features1, features2);

% Retrieve locations of matched points for each image.
matchedPoints1 = valid_points1(indexPairs(:, 1), :);
matchedPoints2 = valid_points2(indexPairs(:, 2), :);
%figure; showMatchedFeatures(img1, img2, matchedPoints1, matchedPoints2, 'montage');
%title('Putative point matches');

%==========================================================================
% Estimate the fundamental matrix (F)
%==========================================================================
% Use the estimateFundamentalMatrix function to compute the F and find the
% inlier points that meet the epipolar constraint. 

[F, inliers] = estimateFundamentalMatrix(...
    matchedPoints1, matchedPoints2,'NumTrials', 4000); 

% Find epipolar inliers. The input 'inlierPts1' represents the coordinates
% of corresponding points in the first view.
inlierPts1 = matchedPoints1(inliers, :);
inlierPts2 = matchedPoints2(inliers, :);
%figure; showMatchedFeatures(img1, img2, inlierPts1, inlierPts2, 'montage');
%title('Epipolar inliers');

%==========================================================================
% Compute epipolar line
%==========================================================================
col = 'bykgrcmw'; % colors for each point-line pair

for i=1:8
    figure(1); title(['Select point #',num2str(i),'/8 with a mouse click'])
    mp = ginput(1)'; % click, then find the conjugate transpose of input point
    plot(mp(1),mp(2),[col(i) 'o'])
    h = labelpoints(mp(1),mp(2),i,'FontSize', 14, 'Color', 'r');
  
    figure(2)
    Drawepipolarline(F,mp,col(i))
end






