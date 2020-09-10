%%-----------------------------------------------------------------------
% Versatile tool to warp PET to template using either the native MRI or
% existing deformation parameters.
%
% Last updated: September 2020
%%-----------------------------------------------------------------------

spm('defaults','PET');

% Prompt to warp via MRI or existing DPs
yn={'Yes','No'};
[warptype,~] = listdlg('PromptString',[{'Do you have the deformation parameters (iy) already?'} {''}],...
                'SelectionMode','multiple',...
                'ListString',yn);

% Prompt for PET regardless of method
pet=spm_select(Inf,'image','Select PET to follow deformations');
sz=size(pet,1); %length defined here

% Procedure for version using DPs
if warptype == 1
    
    % Prompt and prep for this method
    dps=spm_select(Inf,'image','Select deformation parameters (iy_*)');
    dps_spm=cellstr(dps);
    wd=spm_select(1,'dir','Select output directory');
    wd_spm=cellstr(wd);
    pet_spm=cellstr(pet);
    
    % Loop through scans
    for i=1:sz

        % Remove spm_input's ",1" here
        temp_dp=cellstr(dps_spm{i}(1:size(dps_spm{i},2)-2));
        temp_pet=cellstr(pets_spm{i}(1:size(pets_spm{i},2)-2));

        % SPM Deformations batch
        matlabbatch{1}.spm.util.defs.comp{1}.def = temp_dp;
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames = temp_pet;
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = wd_spm;
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {'/home/jagust/ucsf_tau_db/data/templates/icbm152.nii'}; %template used
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.push.prefix = 'w'; %will have this prefix
        spm_jobman('run',matlabbatch);
    
    end
%Procedure for version using MRIs    
elseif warptype == 2
        
    mri=spm_select(Inf,'image','Select MRI (nu.nii) to warp');

    for i=1:sz

        tempmri=mri(i,:);
        temppet=pet(i,:);

        % Define paths for an alternative way to remove spm_input's ",1"
        [p,f,e]=spm_fileparts(tempmri);
        [pp,ff,ee]=spm_fileparts(temppet);

        clear matlabbatch

        % SPM Segment batch to warp MRIs
        % Will output native and template space c1, c2, c3, and DPs
        matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(tempmri);
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,1'};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,2'};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,3'};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,4'};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,5'};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,6'};
        matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.025 0.1];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        
        % SPM Deformations to apply deformation parameters from previous
        % part of batch to the PET
        matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.def = cellstr(strcat(p,'/y_',f,e));
        matlabbatch{2}.spm.util.defs.comp{1}.inv.space = cellstr(strcat(p,'/',f,e));
        matlabbatch{2}.spm.util.defs.out{1}.push.fnames = cellstr(strcat(pp,'/',ff,ee));
        matlabbatch{2}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{2}.spm.util.defs.out{1}.push.savedir.savepwd = {''};
        matlabbatch{2}.spm.util.defs.out{1}.push.fov.file = {'/home/jagust/ucsf_tau_db/data/templates/icbm152.nii'}; %same template
        matlabbatch{2}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{2}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
        matlabbatch{2}.spm.util.defs.out{1}.push.prefix = 'w';
        spm_jobman('run',matlabbatch);

    end
end