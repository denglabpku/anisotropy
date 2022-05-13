function [ LagTime ] = InferFrameRateFromName( FileName )
%INFERFRAMERATEFROMNAME Infer the lag time from the filename
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3

% typical file name:
%   mESC_C59_Halo-mCTCF_16.5Hz_pooled_QC_CD2

% 74 Hz is special:
if ~isempty(strfind(FileName, '_74Hz'))
    LagTime = 0.0135;
elseif ~isempty(strfind(FileName, '_26Hz'))
    LagTime = 0.0385;
elseif ~isempty(strfind(FileName, '_50Hz'))
    LagTime = 0.02;
else
    % Find the location of the frame rate:
    StringLoc = strfind(FileName, 'Hz');

    % So need the 4-5 elements before Hz:
    PotentialFrameRate = FileName(StringLoc-4 : StringLoc-1);

    % now remove any non digits:
    NotFinished = true;
    while NotFinished == true
        if isstrprop(PotentialFrameRate(1), 'digit')
            % it's a number, your job is done
            NotFinished = false;
        else
            % it's not a number, get rid of it:
            PotentialFrameRate(1) = [];
        end
    end
    % convert from string to number
    NumericFrameRate = str2num(PotentialFrameRate);

    % find the lag time in seconds:
    LagTime = 1/NumericFrameRate;
end

end

