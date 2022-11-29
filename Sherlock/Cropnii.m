%Cut the files by running function in sherlock_avg

subjects = {'s1', 's2', 's3', 's4', 's6', 's7', 's9', 's10', 's11', 's12', ...
    's13', 's14', 's15', 's16', 's17'};

% set path to where files are
data_path = 'Users/lucychang/Desktop/Research/sherlock_fmri'; 
cd(data_path);
for s = 1:length(subjects)
    sub = subjects{s};
    filename = [data_path '/' sub '.nii'];
    newname=[sub '_trimmed.nii'];
    sherlock_avg(filename,newname)   
end
