function Node_Import_YZL(FileName,iStudy,varName)
% NODE_IMPORT: Update a brainstorm file with a variable coming from the Matlab base workspace.

% Get node information
iItem    = iStudy;
% Get full filename
[FileName, FileType, isAnatomy] = file_fullpath( FileName );

% Get variable from workspace
[value, varname] = in_matlab_var(varName, 'struct');
if isempty(value)
    return
end

% History: File imported from Matlab variable
if ismember(lower(FileType), {'anatomy', 'scalp', 'outerskull', 'innerskull', 'cortex', 'other', 'channel', 'headmodel', 'data', 'rawdata', 'results', 'kernel', 'pdata', 'presults',  'noisecov', 'dipoles', 'timefreq', 'ptimefreq', 'spectrum', 'matrix'})
    value = bst_history('add', value, 'import', ['Imported from Matlab variable: ' varname]);
end

% Progress bar
bst_progress('start', 'Import from workspace variable', 'Saving file...');
% Save file
save(FileName, '-struct', 'value');
% Reload target subject or study
if isAnatomy
    db_reload_subjects(iItem);
else
    db_reload_studies(iItem);
end
bst_progress('stop');
disp(['BST> File imported from ''' varname '''.']);

% Unload all datasets (safer)
bst_memory('UnloadAll', 'Forced');
% Save database
db_save(); 