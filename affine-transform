%-----------------------------------------------------------------------
% Adapted from /home/jagust/UCSF/rablabtools/AffineWarping.m by Amelia 
% on 6/28/19 to cycle through multiple subjects and MRIs more efficiently,
% especially to transform scans for participants with different numbers of
% scans.
%
% spm SPM - SPM12 (7219)
%-----------------------------------------------------------------------

% grabbing MRIs from all subjects to be transformed
vols_mri = spm_select(Inf,'image', 'Select all MRI images to affine normalize');
vols_mri_spm = cellstr(vols_mri);

waitfor(msgbox({'All the PET scans must' 'be co-registered to the respective' 'MRIs before running this.' 'Have fun!'},'IMPORTANT','warn'));

% loop through MRIs, choose PET to warp for each
for a = 1:size(vols_mri_spm,1)
    
	clear matlabbatch2

    % select PET images
	[p,f,e]=spm_fileparts(vols_mri_spm{a});
	str=sprintf('Select the PET to warp via %s!\n',f);
	vols_pet1 = spm_select(Inf,'image', str); %% images 
	vols_pet1_spm=cellstr(vols_pet1); %% same

	b = size(vols_pet1_spm,1);

    % prep all subject-specific images for the batch
	temp_mri=cellstr(vols_mri_spm{a});
	imgs=vertcat(temp_mri,vols_pet1_spm);

	spm('defaults','PET');

    % set up and run normalization batch through SPM
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.subj(1).source = temp_mri; %uses MRI to define parameters
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.subj(1).wtsrc = '';
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.subj(1).resample = imgs; %transforms MRI and all PET
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.template = {'/srv/local/matlab-tools/spm/spm12/toolbox/OldNorm/T1.nii,1'};
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 0;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-100 -130 -80
		                                                 100 100 110];
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.vox = [1 1 1];
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1; %trilinear interpolation
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
	matlabbatch2{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w_affine_'; %defines new filenames
	spm_jobman('run',matlabbatch2);



end

clear;
