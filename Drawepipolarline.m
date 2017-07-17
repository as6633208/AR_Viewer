% This function computes the eipolar line in second view from fundamental
% matrix (F) and point m1 in frist view. The line is drawn using color col.
%
%    Input:
%       - F, a 3-by-3 fundamental matrix, if m1 represent a point in the
%       first iamge I1 that corresponds to m2, a point in the second image
%       I2, then: m2'*F*m1=0
%       
%       - point, an M-by-2 matrix, where each row contains the x and y
%       corredinates of a point in image. M represents the number of
%       points.
%
%    Output:
%       - lines, an M-by-3 matrix, where each row must be in the format
%       [a,b,c]. M represents the number of lines.

function Drawepipolarline(F,m1,col)
hold on
if length(m1) == 2   
    m1 = [m1; 1];
end

epipLine = F*m1;    % computes line
ax  = axis;         % returns the x-axis and y-axis limits for the current axes. 
x   = ax(1:2)';     % Find the conjugate transpose of ax (extract the first through the second elements).   
a   = epipLine(1);  % A line in 2D is represented by the homogeneous 3-vector
b   = epipLine(2); 
c   = epipLine(3); 
y   = -(c+a*x)/b; 

plot(x,y,x,y,col)





