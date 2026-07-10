function [sldProfiles, resampledLayers, pIntervals] = alignALProfiles(geometry, modelType, sldProfiles, resampledLayers, pIntervals)
% Aligns the A/L SLD profiles so that the substrates line up by padding the
% start of any shorter than the longest profile. Also adjusts resampled
% layers and prediction interval slds as necessary.
%
% Parameters
% ----------
% geometry : geometryOptions
%    The geometry.
% modelType : modelTypes
%    The model type. 
% sldProfiles : cell
%    The sld profiles.
% resampledLayers : cell
%    The resampled layers. 
% pIntervals : cell
%    The slds in the prediction intervals.
%
% Returns
% -------
% sldProfiles : cell
%    The sld profiles, adjusted if necessary.
% resampledLayers : cell
%    The resampled layers, adjusted if necessary. 
% pIntervals : cell
%    The slds in the prediction intervals, adjusted if necessary.

if ~strcmpi(geometry,'air/substrate') || strcmpi(modelType,'custom xy')
    return
end

% Find the length of the longest profile.
lengths = zeros(size(sldProfiles));
for i = 1:numel(sldProfiles)
    lengths(i) = size(sldProfiles{i}, 1);
end

% Get max length and its index
[maxLen, maxPos] = max(lengths, [], "all");
maxXValue = sldProfiles{maxPos}(end,1);

% Get the longest profile...
maxX = sldProfiles{maxPos}(:,1);

% Pad the start of any profiles that are shorter than this
for i = 1:numel(sldProfiles)
    thisSLD = sldProfiles{i};
    thisLen = size(thisSLD,1);
    if thisLen < maxLen
        diffLen = maxLen - thisLen;
        pad = zeros(diffLen,1);
        newY = [pad ; thisSLD(:,2)];
        sldProfiles{i} = [maxX(:,1) newY(:)];
        
        % For resampled layers, the pad is just one big layer at the start
        thisResam = resampledLayers{i};
        if any(thisResam,'all')      % not all zeros
            totLength = sum(thisResam(:,1));
            padLength = maxXValue - totLength;
            resamPad = [padLength 0 0];
            resampledLayers{i} = [resamPad ; thisResam]; 
        end

        if ~isempty(pIntervals)
            pIntervals{i} = [zeros(5, diffLen), pIntervals{i}];
        end
    end
end
end
