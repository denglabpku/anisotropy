function [ Finish ] = EditRunInputFile_for_batch( FileName, inputfile, outputfile, LagTime )
%EDITRUNINPUTFILE_FOR_BATCH Edit the RunInputFile to run vbSPT in batch
%mode
%   written by Anders Sejr Hansen (AndersSejrHansen@post.harvard.edu;
%   @Anders_S_Hansen; https://anderssejrhansen.wordpress.com)
%   License: GNU GPL v3

%FileName = 'vbSPT_RunInputFileBatch.m';

% Read txt into cell called LineCell, do this line-by-line
LineCell = regexp( fileread(FileName), '\n', 'split');

% the line numbers to change are:
%   inputfile is line number 9
%   outputfile is line number 18
%   timestep is line number 22

NewLineCell = LineCell;

% change the key lines:
NewLineCell{9} = ['inputfile = ', '''', inputfile, ''';'];
NewLineCell{18} = ['outputfile = ', '''',outputfile, ''';'];
NewLineCell{22} = ['timestep = ', num2str(LagTime), '; %timesep in [s]']; 


% now overwrite the existing file:
fid = fopen(FileName, 'w');
fprintf(fid, '%s\n', NewLineCell{:});
fclose(fid);

% OK, all done
Finish = 1;

end

