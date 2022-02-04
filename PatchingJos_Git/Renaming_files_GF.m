% Get all abf files in the current folder
files = dir('*.abf');
% A = datetime(files.date, 'InputFormat','yyyy-MM-dd HH:mm:ss');

files = files(~[files.isdir]);
[~,idx] = sort([files.datenum]);
files = files(idx);

% Loop through each
for id = 1:length(files)
    filename = files(id).name;
    newfilename = sprintf([filename(1:end-4), '_GapFree', '.abf'],id);
    movefile(filename,newfilename);
end
