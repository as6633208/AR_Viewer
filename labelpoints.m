function [h, ext] = labelpoints (xpos, ypos, labels, varargin)

% Check Class of 'labels'
    %If 'labels' are numerical, convert to cell
    if isnumeric(labels) == 1
        labels = num2cell(labels); 
    end

    % if 'labels' are char, convert to cell
    if ischar(labels)
        labels = cellstr(labels);
    end
    
% if all labels share the same xpos or ypos (only 1 value entered in 1 of the position vectors)
    if length(xpos)==1 && length(ypos)>1
        xpos = repmat(xpos, size(ypos));
    elseif length(ypos)==1 && length(xpos)>1
        ypos = repmat(ypos, size(xpos));
    end
    
% if only one label is entered for all points, replicate it
    if length(labels)==1 && length(xpos) > 1
        labels = repmat(labels, [1, length(xpos)]);
    end
    
% ensures xpos, ypos, and labels are all row vectors 
    if iscolumn(xpos);      xpos = xpos';       end
    if iscolumn(ypos);      ypos = ypos';       end
    if iscolumn(labels);    labels = labels';   end

%By this point, xpos, ypos and labels should all be same length EXCEPT if optional input 'stacked' is used
% indicate error otherwise.
    if isequal(length(xpos), length(ypos), length(labels)) == 0 && sum(strcmp('stacked', varargin))==0
        error('xpos, ypos, and labels must all be the same length unless using one input for labels.')
    end
    
%if an 'inf' value is entered, this will eliminate that entire point and label
    xinf = find(xpos==inf);
    yinf = find(ypos==inf);
    findinf = [xinf yinf];
    if ~isempty(findinf)
        xpos(findinf)=[];
        ypos(findinf)=[];
        labels(findinf) = [];
    end       
    
%Validate inputs and optional parameters
  %for help, see https://www.mathworks.com/help/matlab/ref/inputparser.addparameter.html#inputarg_validationFcn
  % and https://www.mathworks.com/help/matlab/ref/validateattributes.html
    validPositions = {'N' 'NE' 'E' 'SE' 'S' 'SW' 'W' 'NW' 'center' 'C'};
    checkPosition = @(x) any(validatestring(x, validPositions));

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'xpos', @isnumeric);
    addRequired(p, 'ypos', @isnumeric);
    addRequired(p, 'labels');
    addOptional(p, 'position', 'NW', checkPosition);
    addOptional(p, 'buffer', 0, @isnumeric);
    addOptional(p, 'adjust_axes', 0, @isnumeric);

    addParameter(p, 'outliers_SD', 3, @isnumeric);
    addParameter(p, 'outliers_Q', 1.5, @isnumeric);
    addParameter(p, 'outliers_N', 1, @isnumeric);
    addParameter(p, 'outliers_lim', [0,0;0,0]);
    addParameter(p, 'outliers_lin', {'SD',1});
    
    addParameter(p, 'stacked', 'down');

    addParameter(p, 'FontSize', 10, @isnumeric);
    addParameter(p, 'Color', 'k');
    addParameter(p, 'rotation', 0, @isnumeric);
    parse(p, xpos, ypos, labels, varargin{:})

%calculate buffer
    a = axis/10;% I've somewhat arbitrarily divided by 10 to make 'buffer' more sensitive
    u1 = 0;     %x offset
    u2 = 0;     %y offset
    
%assign position
    switch upper(p.Results.position) 
        case 'E',       va = 'middle'; ha = 'left';         u1 = a(2)-a(1);         
        case 'W',       va = 'middle'; ha = 'right';        u1 = (a(2)-a(1))*-1;
        case 'N',       va = 'bottom'; ha = 'center';       u2 = a(4)-a(3);
        case 'S',       va = 'top';    ha = 'center';       u2 = (a(4)-a(3))*-1;
        case 'NE',      va = 'bottom'; ha = 'left';         u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))/2;
        case 'NW',      va = 'bottom'; ha = 'right';        u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))/2;
        case 'SE',      va = 'top';    ha = 'left';         u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))*-0.5;
        case 'SW',      va = 'top';    ha = 'right';        u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))*-0.5;
        case {'CENTER', 'C'},  va = 'middle'; ha = 'center'; 
    end
    
%Factor in buffer
    u1 = u1*p.Results.buffer;   %applied to X
    u2 = u2*p.Results.buffer;   %applied to Y
    
%adjust u1, u2 if rotation is enabled 
%  rotation centers text on data point no matter what the position input is. 
%  so here we increase the distance from the point a tad to compensate for the lack of true alignment. 
    if p.Results.rotation~=0           
        factor = 1.2; %arbitrary, 1=no change, <1 means shift toward dot, >1 means shift away from dot. 
        u1 = u1*factor;
        u2 = u2*factor;
        %if we are rotating the text, we'll initially plot it centered - this must happen before plotting labells.
        va = 'middle'; ha = 'center';
    end
    
% Index (1/0) of all paired data (ie, if 1 or both coordinates are NaN, pairIDX(i) is 0.
% This may be used in identifying outliers but should NOT be used to calculate stats on xpos or ypos
    pairIDX = ~isnan(xpos) & ~isnan(ypos);

%% stacked text section (should come before outlier section)
if sum(strcmp('stacked', varargin))==1
    %if optional input 'stacked' is being used, make sure user isn't confused by using option incompatable params
    %if more than 1 xpos or ypos is entered...
    if length(xpos)>1 || length(ypos)>1                     
        warning('Only the first xpos and ypos will be used to initiate stacked text');
    end
    %if outliers are entered in input, remove them from varargin
    if ~isempty(cell2mat(regexp(varargin(cellfun(@ischar, varargin)), 'outlier')))      
        warning('Cannot use stacked and outlier input parameters at the same time.  Outliers input will be ignored');
        tempvarargin = varargin;
        tempvarargin(cellfun(@isnumeric, tempvarargin)) = {'temp_replace'};    %numeric values prevent use of regexp()
        varargIDX = find(~cellfun(@isempty,regexp(tempvarargin, 'outlier')));
        varargin([varargIDX, varargIDX+1]) = [];    %removes outlier input and its associated param from vararin
    end
    
    %internal parameter defaults
    spacing = 1;
    stacktype = lower(p.Results.stacked); 
    
    %detect if user added optional spacing parameter (example: 'down_1.5')
    %if detected, this separates the stacktype from spacing 
    if ~isempty(strfind(stacktype, '_'))
        spacing = str2double(stacktype(strfind(stacktype, '_')+1:end));
        stacktype = stacktype(1:strfind(stacktype, '_')-1);
    end
    
    %Check that user entered a valid stacktype
    if ~any(strcmp(stacktype, {'up', 'down', 'left', 'right'}))
        error('Text stacking options are:  up, down, left, or right (not case sensitive)');
    end        
    
    %clear xpos and ypos vectors after initial points
    nlabs = length(labels); %number of labels
    if nlabs>1              %this is needed if user only has 1 element in stack. 
        xpos(min(2, nlabs):nlabs) = nan;
        ypos(min(2, nlabs):nlabs) = nan;
    end
    
   %get xpos and ypos for all other labels
    for s = 2:nlabs
        
        %Temperarily plot the s-1 label, get its position, then delete it.
        temphand = text(xpos(s-1)+u1 , ypos(s-1)+u2, labels(s-1), 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize);
        tempextnt = get(temphand, 'extent');
        delete(temphand)
        
        %Calculate xpos and ypos for label S (this came from stackedtext.m which is now obsolete)
        switch stacktype
            case 'down'
                xpos(s) = xpos(s-1);
                ypos(s) = ypos(s-1) - tempextnt(4)*spacing;
            case 'right'
                ypos(s) = ypos(s-1);
                xpos(s) = xpos(s-1) + tempextnt(3)*spacing;
            case 'up'
                xpos(s) = xpos(s-1);
                ypos(s) = ypos(s-1) + tempextnt(4)*spacing;
            case 'left'
                ypos(s) = ypos(s-1);
                xpos(s) = xpos(s-1) - tempextnt(3)*spacing;
        end %switch  
        
    end %for s = 2
    
end %if sum()
    

%% Outlier section
%If outliers parameters are selected
    outlierNames = {'outliers_SD', 'outliers_Q', 'outliers_N', 'outliers_lim', 'outliers_lin'}; %add new outlier options here
    outlier_flag = false;
    %identify if/which outlier inputs is (or isn't) present
    for c = 1:length(outlierNames)
        if sum(strcmp(varargin, outlierNames{c}))>0;
            outlier_flag = true;  
            outliertype = varargin{strcmp(varargin, outlierNames{c})};  %cell naming outlier type
        end
    end
    %executes chosen outlier type.  The idea is that each case (or type) serves the purpose of identifying what outliers to keep. 
    %so the output of each subsection is the 'outlier_idx' variable which is an index of all labels'
    if outlier_flag 
        switch outliertype
            
            case 'outliers_SD'
                SDs = p.Results.outliers_SD(1);
                % if user specified the 'center' of her data:
                if length(p.Results.outliers_SD) > 1
                    Xcnt = p.Results.outliers_SD(2);
                    Ycnt = p.Results.outliers_SD(3);
                else % if user did not specify center of her data:
                    Xcnt = nanmean(xpos);
                    Ycnt = nanmean(ypos);
                end
                outlier_idx = logical(abs(xpos - Xcnt) > SDs*nanstd(xpos)  |  abs(ypos - Ycnt) > SDs*nanstd(ypos)); %index of outliers
            
            case 'outliers_Q'
                xbounds = [prctile(xpos,25) - p.Results.outliers_Q * iqr(xpos) , prctile(xpos, 75) + p.Results.outliers_Q * iqr(xpos)];   %[lower upper] bounds of outliers
                ybounds = [prctile(ypos,25) - p.Results.outliers_Q * iqr(ypos) , prctile(ypos, 75) + p.Results.outliers_Q * iqr(ypos)];   %[lower upper] bounds of outliers
                outlier_idx = logical(ypos<ybounds(1) | ypos>ybounds(2) |  xpos<xbounds(1) | xpos>xbounds(2));
            
            case 'outliers_lim'
                %assign limits and qualifier
                limvars = p.Results.outliers_lim; %ie, {[5 5; 5 5]; 'or'}
                if iscell(limvars)
                    lims = limvars{1}; 
                    if size(limvars,1) == 1 %default
                        qualifier = 'or';
                    else  qualifier = lower(limvars{2});
                    end
                else
                    lims = limvars;
                    qualifier = 'or';
                end

                if size(lims,1) == 1
                    lims = [lims;lims];
                end

                x_outliers = xpos<lims(1,1) | xpos>lims(1,2);  %idx of x points outside of safe zone
                y_outliers = ypos<lims(2,1) | ypos>lims(2,2);  %idx of y points outside of safe zone
                switch qualifier
                    case 'or'
                        outlier_idx = x_outliers | y_outliers;
                    case 'and'
                        outlier_idx = x_outliers & y_outliers;
                    case 'invert'
                        x_outliers = xpos>lims(1,1) & xpos<lims(1,2);
                        y_outliers = ypos>lims(2,1) & ypos<lims(2,2);
                        outlier_idx = x_outliers & y_outliers;
                    otherwise
                        error('The inputs you entered for Outliers_lim wasn''t recognized.')
                end
                
            case 'outliers_N'
                Npts = p.Results.outliers_N(1);
                % if user specified the 'center' of her data:
                if length(p.Results.outliers_N) > 1
                    Xcnt = p.Results.outliers_N(2);
                    Ycnt = p.Results.outliers_N(3);
                else % if user did not specify center of her data:
                    Xcnt = nanmean(xpos);
                    Ycnt = nanmean(ypos);
                end
                if p.Results.outliers_N<1;  
                    N = ceil(length(xpos(pairIDX)) * Npts); 
                else
                    N = min(Npts, length(xpos(pairIDX)));        %ensures that user cannot label more outliers than coordinates.
                end
                 meanpoint = repmat([Xcnt Ycnt], [length(xpos),1]);
                 paired = horzcat(xpos', ypos');
                 distances = (((meanpoint(:,1)-paired(:,1)).^2)  +  ((meanpoint(:,2)-paired(:,2)).^2)).^(1/2);       %all distances from mean
                 [sorted, idx] = sort(distances, 'descend');
                 idx = idx(~isnan(sorted)); %this is to ignore any NaN values in xpos or ypos that would cause an nan distance. 
                 outlier_idx = false(1,length(xpos));
                 outlier_idx(idx(1:N))=1; 
                 
            case 'outliers_lin'
                %user can either enter {slope, yint, outlier type, threshold} -or- {outlier type, threshold}
                %here we control for what user entered or didn't enter
                linvars = p.Results.outliers_lin; %ie, {1,1.1,'sd',1} or {'sd',1}
                if isnumeric(linvars{1}) && ~isempty(linvars{1}) %if user entered own slope and y-int
                    slope = linvars{1};
                    yint = linvars{2};
                    outtype = upper(linvars{3}); %outlier type (sd or n)
                    outthresh = linvars{4};      %threshold (values > threshold are outliers)
                else
                    %calculated slope and y-int of the (x,y) data.
                    if isempty(linvars{1}); linvars = linvars([3,4]); end
                    slope = nansum((xpos-nanmean(xpos)).*(ypos-nanmean(ypos))) / nansum((xpos-nanmean(xpos)).^2);       
                    yint = nanmean(ypos) - (slope * nanmean(xpos));                     %refline(slope,yint)   %To test
                    outtype = upper(linvars{1});
                    outthresh = linvars{2};
                end
                %now calculate residuals (linear estimate - y val)^2
                Yest = slope*xpos + yint;
                resid = (Yest - ypos).^2;
                %now sort residuals > to < similar to 'outliers_N'
                [sorted, idx] = sort(resid, 'descend');
                idx = idx(~isnan(sorted)); %this is to ignore any NaN values in xpos or ypos that would cause an nan distance. 
                %finally, chose the outliers based on chosen method
                if strcmp(outtype, 'SD')
                    outlier_idx = idx(1:sum(nanstd(sorted)*outthresh <= sorted));     %figure; plot(sorted, 'o'); addhorz(nanstd(sorted)*outthresh);  %to view sorted outliers and treshold
                elseif strcmp(outtype, 'N')
                    if outthresh<1;  
                        N = ceil(length(idx) * outthresh); 
                    else
                        N = min(outthresh, length(idx));        %ensures that user cannot label more outliers than coordinates.
                    end
                    outlier_idx = idx(1:N);
                end
        end %outlier switch
        
        xpos = xpos(outlier_idx);               
        ypos = ypos(outlier_idx);
        labels = labels(outlier_idx);

%         if any(outlier_idx) == 0;           %dispay msg if there are no outliers to label
%             mfile = [mfilename,'.m'];
%             disp(['There are no outliers to label in ', mfile,'.'])
%             disp('Change outlier value for less sensitivity; See help file.'); 
%         end
    end %outlier_flag
%%

%Label points
    hand = text(xpos+u1 , ypos+u2, labels, 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize, 'color', p.Results.Color);
    extnt = get(hand, 'extent');
    
%Rotate text if specified
    if p.Results.rotation~=0           %if rotation input is something other than 0  (use to be:  sum(strcmp(varargin, 'rotation')) == 1 )
        xl = xlim;      yl = ylim;                          %In case text rotation auto adjusts axes.
        curr_extent = get(hand, 'extent');                     %Need to store current center point of all labels since text rotation relocates position
        if iscell(curr_extent); curr_extent = cell2mat(curr_extent); end
        hold on
        curr_position = [curr_extent(:,1)+(curr_extent(:,3)/2),curr_extent(:,2)+(curr_extent(:,4)/2)];          %uses extent to locate center of label
        set(hand, 'rotation', p.Results.rotation, 'VerticalAlignment','middle', 'HorizontalAlignment','center');  	%note: text rotation changes alignment which is why they need to be centered back to specifications.
        for i = 1:length(hand)                                 %after rotation, reposition labels back to desired location 
            set(hand(i), 'position', curr_position(i,:))
        end
        set(gca, 'xlim', xl); set(gca, 'ylim', yl);         %In case text rotation auto adjusts axes.
    end     
    
%Determine if any labels go beyond axis limits and adjust if desired  (adjust_axes = 0 or 1)
    if p.Results.adjust_axes == 1   &&   ~isempty(hand)    
        x_adj = sign(u1+0.0000001);                 %the addition is to avoid '0'
        y_adj = sign(u2+0.0000001);                 %the addition is to avoid '0'

        labelextent = get(hand, 'extent');
        if isequal(class(labelextent),'cell')
           labelextent = cat(1, labelextent{:});
        end
        xl = xlim;      yl = ylim;
        lablimX = [min(labelextent(:,1)), max(labelextent(:,1)+(labelextent(:,3).*x_adj))] +u1;
        lablimY = [min(labelextent(:,2)), max(labelextent(:,2)+(labelextent(:,4).*y_adj))] +u2;

        xlim([min(min(xl), min(lablimX)), max(max(xl), max(lablimX))])
        ylim([min(min(yl), min(lablimY)), max(max(yl), max(lablimY))])
    end
    
%Turn off Latex interpreter to avoid subscripts with an underscore is used in a label
    set(hand, 'interpreter', 'none')
    
%Outputs
if nargout > 0 
    h   = hand;
    ext = extnt;
end
     
end