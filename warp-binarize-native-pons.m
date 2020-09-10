%%---------------------------------------------------------------------------
% Project-specific script. We had FDG images processed with an MRI-free
% pipeline that we wanted to process using a procedure that was as close as
% possible to our in-house FDG processing.
%
% This script was written to warp a specific template pons ROI to native space
% (Deformations) then binarize (ImCalc), using the SPM native c2 (WM
% probability) as an extra security measure. We can then use this native space
% pons ROI to renormalize our FDG images for analysis.
%
% Written by Amelia March 2020
%%-----------------------------------------------------------------------------

spm('defaults','PET');

% Collect all necessary images
dps=spm_select(Inf,'image','Select deformation parameters (y_*)');
nus=spm_select(Inf,'image','Select all nu.nii images');
c2s=spm_select(Inf,'image','Select all c2 images');
wd=spm_select(1,'dir','Select output directory for non-c2-masked ponses');
od=spm_select(1,'dir','Select output directory for final ponses');

% Ensure first outputted images are located in the correct directory
cd(wd);

[sz,~]=size(dps); %used as default length for all loops

% Prep images
dps_spm=cellstr(dps);
nus_spm=cellstr(nus);
c2s_spm=cellstr(c2s);

% Use regex to find 4-digit IDs in each DP filename
ids=regexp(dps_spm,'\d{4}','match','once'); %creates nx1 cell matrix of 4-digit IDs

tpons={'/home/jagust/astrom/pcc_project/ADNI/scripts/rpons/mean_w_NN_rpons_allbin_bin0p1.nii'}; %template pons

% Loop through each subject
for i=1:sz

    clear matlabbatch
    clear matlabbatch2

    % These SPM batches require specific formats (i.e. excluding the ",1" included at the end of all spm_input strings)
    dptemp=cellstr(dps_spm{i}(1:size(dps_spm{i},2)-2));
    nutemp=cellstr(nus_spm{i}(1:size(nus_spm{i},2)-2));
    c2temp=cellstr(c2s_spm{i}(1:size(c2s_spm{i},2)-2));
    idtemp=ids{i};
    prefix=strcat('w_',idtemp,'_'); %define prefix added once ponses are warped

    % SPM Deformations batch to warp template pons to native space
    matlabbatch{1}.spm.util.defs.comp{1}.def = dptemp; %dp
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames = tpons; %template pons
    matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.savepwd = 1;
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = nutemp; %native MRI as FOV
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix = prefix;
    spm_jobman('run',matlabbatch);

    output=strcat(idtemp,'_rpons_mask_c2'); %define final binarized pons
    ppons=cellstr(strcat(wd,'/',prefix,'mean_w_NN_rpons_allbin_bin0p1.nii')); %define native pons just created
    temp2=vertcat(ppons,c2temp); %combine native pons and c2 for batch

    % SPM ImCalc batch to binarize native pons (at pons>0.5) using c2 also (at c2>0.2)
    matlabbatch2{1}.spm.util.imcalc.input = temp2; %pons and c2
    matlabbatch2{1}.spm.util.imcalc.output = output;
    matlabbatch2{1}.spm.util.imcalc.outdir = cellstr(od);
    matlabbatch2{1}.spm.util.imcalc.expression = '(i1>0.5).*(i2>0.2)';
    matlabbatch2{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch2{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch2{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch2{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch2{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch2);

end

clear;
