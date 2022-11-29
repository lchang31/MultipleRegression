function sherlock_avg(filename,newname)

[fpath filename_noext] = fileparts(filename);
nii = load_untouch_nii(filename);
dsize = size(nii.img);
nii.img = avg_img(nii.img);
nii.hdr.dime.dim(5) = size(nii.img,4);
new_filename = newname;
save_untouch_nii(nii,fullfile(fpath,new_filename));
% return the number of TR
niiFrame = get_nii_frame(newname);
X = ['number of TR = ', num2str(niiFrame)];
disp(X);

% takes the average of TR
function data = avg_img(data)
    i = 1;
    k = 1;
    while i < 1976
        data(:,:,:,k) = (data(:,:,:,i) + data(:,:,:,i+1))/2;
        i = i+2;
        k = k+1; 
    end
    data(:,:,:,988+1:end) = [];
end    
end
