% Get all abf files in the current folder
files = dir('*.abf');
% A = datetime(files.date, 'InputFormat','yyyy-MM-dd HH:mm:ss');

files = files(~[files.isdir]);
[~,idx] = sort([files.datenum]);
files = files(idx);

% Loop through each
for id = 1:length(files)
    filename = files(id).name;
%     newfilename = sprintf([filename(1:end-4), '_GapFree', '.abf'],id);
  newfilename = sprintf([filename(1:4), '_', filename(5:6), '_', filename(7:end)],id);
% newfilename = sprintf([filename(1:12), '03', filename(15:end-6), sprintf('%02d', id), '.abf'],id);
%   newfilename = sprintf([filename(1:13), '0', filename(15:end)],id);
     movefile(filename,newfilename);
end
