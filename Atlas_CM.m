%% Apply atlas based parcellation on MNI fuctional data

% func_path = full path to MNI 4D preprocessed fMRI image, filtered_func_data_clean_MNI
% atlas_path = full path to atlas image
% outdir = full path to output directory


function Atlas_CM(func_path,subjdirs,atlas_path,outdir)
    cd(func_path);
    atlas=niftiread(atlas_path);
    atlas_cortex_vector=atlas(:);
    num_Regions_cortex=length(unique(atlas_cortex_vector(atlas_cortex_vector~=0)));
    rois=unique(atlas_cortex_vector(atlas_cortex_vector~=0));
    for j=1:length(subjdirs)
        sub=subjdirs{j};
        disp(sub)
        if isfile([outdir '/' sub '.mat'])==0
            func_data=niftiread([sub '/MNINonLinear/Results/rsfMRI_AP/rsfMRI_AP_hp2000_clean.nii.gz']);
            num_Tpoints=size(func_data,4);
            func_data_vector=reshape(func_data,[],num_Tpoints);
            mean_ts_cortex=zeros(num_Regions_cortex,num_Tpoints);
            for k=1:num_Regions_cortex
                roi=rois(k);
                vox=find(atlas_cortex_vector==roi);
                tmpData=func_data_vector(vox,:);
                mean_ts_cortex(k,:)=mean(tmpData);
            end
            CM_cortex=corr(mean_ts_cortex');
            z_CM_cortex= .5*(log(1+CM_cortex) - log(1-CM_cortex));
            atlas_z_file_path=[outdir '/' sub '.mat'];
            save(atlas_z_file_path,'z_CM_cortex');
        end
    end
    
end