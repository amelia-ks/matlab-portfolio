%%-----------------------------------------------------------------------
% Written by Amelia to automate longitudinal PET processing using an eroded
% white matter mask as the reference region.
%   1. MRIs are run through pairwise longitudinal registration in SPM to
%      create a "midpoint" MRI
%   2. Midpoint MRI is segmented in SPM to create c2 probability image
%   3. PET are co-registered to midpoint MRI to ensure proper extraction
%   4. c2 is binarized and eroded to create reference region
%   5. PET images are renormalized by mean values in reference region
%
% This script assumes that:
%  1. All PET are coregistered to their respective MRI
%  2. PET images have different filenames for their respective timepoints
%%-----------------------------------------------------------------------


% Adding path to ErodeMask.m script just in case
addpath('/home/jagust/UCSF/pet_by_hand/longitudinal/scripts');

% Ask for output directory
wd=spm_select(1,'dir','Select output directory');
wd_spm=cellstr(wd);

cd (wd);

% Get MRIs and dates to calculate time between scans
vols_mri=spm_select(2,'image','Select the 2 MRIs in order');
vols_mri_spm=cellstr(vols_mri);

waitfor(msgbox({'The following dates must' 'be in the format MM/dd/yyyy' 'Have fun!'}));

mri_date1=spm_input('Date of MRI 1?',1,'s');
mri_date2=spm_input('Date of MRI 2?',2,'s');

% Find out what PET we need to process with the MRI
str={'FTP only','PIB only','FTP and PIB','No PET, just MRI'};
                 [pet_mods,v] = listdlg('PromptString',[{'What PET do you need to run longitudinally?'} {''}],...
                'SelectionMode','single',...
                'ListString',str);

% Ask for PET scans according to previous answer
if pet_mods==1||3
    
    vols_ftp=spm_select(2,'image','Select the 2 FTP scans in order');
    vols_ftp_spm=cellstr(vols_ftp);
    
end

if pet_mods==2||3

	vols_pib=spm_select(2,'image','Select the 2 PIB scans in order');
	vols_pib_spm=cellstr(vols_pib);

end

% Calculates number of years between dates for Jacobian creation
numdays=datenum(mri_date2)-datenum(mri_date1);
numyears=numdays/365.25;

spm('defaults','PET');

% Run longitudinal registration to get midpoint MRI and Jacobians
matlabbatch{1}.spm.tools.longit.pairwise.vols1 = vols_mri_spm(1);
matlabbatch{1}.spm.tools.longit.pairwise.vols2 = vols_mri_spm(2);
matlabbatch{1}.spm.tools.longit.pairwise.tdif = numyears;
matlabbatch{1}.spm.tools.longit.pairwise.noise = NaN;
matlabbatch{1}.spm.tools.longit.pairwise.wparam = [0 0 100 25 100];
matlabbatch{1}.spm.tools.longit.pairwise.bparam = 1000000;
matlabbatch{1}.spm.tools.longit.pairwise.write_avg = 1;
matlabbatch{1}.spm.tools.longit.pairwise.write_jac = 1;
matlabbatch{1}.spm.tools.longit.pairwise.write_div = 1;
matlabbatch{1}.spm.tools.longit.pairwise.write_def = 1;
spm_jobman('run',matlabbatch);

[~,f,e]=spm_fileparts(vols_mri(1,:));

% Define midpoint MRI
avg_mri=strcat(wd,'/avg_',f,e);
avg_mri_spm=cellstr(avg_mri);

% Define jacobians for later

clear matlabbatch

% Segment midpoint MRI to get c1, c2, and c3
matlabbatch{1}.spm.spatial.preproc.channel.vols = avg_mri_spm;
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/usr/local/matlab-tools/spm/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
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
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm_jobman('run',matlabbatch);

% Define c1, c2 and c3 images
[p,f,e]=spm_fileparts(avg_mri);

avg_c1=strcat(p,'/c1',f,e);
avg_c1_spm=cellstr(avg_c1);
avg_c2=strcat(p,'/c2',f,e);
avg_c2_spm=cellstr(avg_c2);
avg_c3=strcat(p,'/c3',f,e);
avg_c3_spm=cellstr(avg_c3);

% Define c2 mask image for later
c2bin=strcat('c2',f,'_bin',e);

% Coregister all PET to midpoint MRI
if pet_mods==1||3

    clear matlabbatch
    
	matlabbatch{1}.spm.spatial.coreg.estwrite.ref = avg_mri_spm;
	matlabbatch{1}.spm.spatial.coreg.estwrite.source = vols_ftp_spm(1);
	matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
	matlabbatch{2}.spm.spatial.coreg.estwrite.ref = avg_mri_spm;
	matlabbatch{2}.spm.spatial.coreg.estwrite.source = vols_ftp_spm(2);
	matlabbatch{2}.spm.spatial.coreg.estwrite.other = {''};
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
	spm_jobman('run',matlabbatch);

	% Define coregistered images
	[~,f,e]=spm_fileparts(vols_ftp_spm{1});
	reg_ftp1=strcat(wd,'/r',f,e);
	reg_ftp1_spm=cellstr(reg_ftp1);
    
    %Name renormalized images for later
    renorm_ftp1=strcat('r',f,'_renormWM',e);

	[p,f,e]=spm_fileparts(vols_ftp_spm{2});
	reg_ftp2=strcat(wd,'/r',f,e);
	reg_ftp2_spm=cellstr(reg_ftp2);
    
    %Name renormalized images for later
    renorm_ftp2=strcat('r',f,'_renormWM',e);
    
end

if pet_mods==2||3

    clear matlabbatch
    
	matlabbatch{1}.spm.spatial.coreg.estwrite.ref = avg_mri_spm;
	matlabbatch{1}.spm.spatial.coreg.estwrite.source = vols_pib_spm(1);
	matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
	matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
	matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
	matlabbatch{2}.spm.spatial.coreg.estwrite.ref = avg_mri_spm;
	matlabbatch{2}.spm.spatial.coreg.estwrite.source = vols_pib_spm(2);
	matlabbatch{2}.spm.spatial.coreg.estwrite.other = {''};
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
	matlabbatch{2}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.interp = 4;
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.mask = 0;
	matlabbatch{2}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
	spm_jobman('run',matlabbatch);

	% Define coregistered images
	[~,f,e]=spm_fileparts(vols_pib_spm{1});
	reg_pib1=strcat(wd,'/r',f,e);
	reg_pib1_spm=cellstr(reg_pib1);
    
    renorm_pib1=strcat('r',f,'_renormWM',e);

	[~,f,e]=spm_fileparts(vols_pib_spm{2});
	reg_pib2=strcat(wd,'/r',f,e);
	reg_pib2_spm=cellstr(reg_pib2);
    
    renorm_pib2=strcat('r',f,'_renormWM',e);
    
end

clear matlabbatch

% Binarize c2 image at 0.01
matlabbatch{1}.spm.util.imcalc.input = avg_c2_spm;
matlabbatch{1}.spm.util.imcalc.output = c2bin;
matlabbatch{1}.spm.util.imcalc.outdir = wd_spm;
matlabbatch{1}.spm.util.imcalc.expression = 'i1>0.01';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);

% Use existing script to erode c2 mask
ErodeMask(c2bin,3);
firsterode=strcat('erode3_',c2bin);
ErodeMask(firsterode,3);

% Remove masks that are now defunct
delete (c2bin);
delete (firsterode);

% Reshape c2 mask for getting mean intensity
eroded_c2=strcat(wd,'/erode3_',firsterode);
Veroded_c2=spm_vol(eroded_c2);
Veroded_c2_read=spm_read_vols(Veroded_c2);
[sz1,sz2,sz3]=size(Veroded_c2_read);
reroded_c2=reshape(Veroded_c2_read,sz1*sz2*sz3,1);

if pet_mods==1||3

	% Reshape PET images for getting mean intensity
	FTP1=spm_vol(reg_ftp1);
	FTP1_read=spm_read_vols(FTP1);
	rFTP1=reshape(FTP1_read,sz1*sz2*sz3,1);

	FTP2=spm_vol(reg_ftp2);
	FTP2_read=spm_read_vols(FTP2);
	rFTP2=reshape(FTP2_read,sz1*sz2*sz3,1);

	% Get mean intensities
	ind_ftp1=find(reroded_c2>0.5 & rFTP1>0);
	meanWM_ftp1=mean(rFTP1(ind_ftp1));
    mWM_ftp1=num2str(meanWM_ftp1);

	ind_ftp2=find(reroded_c2>0.5 & rFTP2>0);
	meanWM_ftp2=mean(rFTP2(ind_ftp2));
    mWM_ftp2=num2str(meanWM_ftp2);

	% Define ImCalc expressions for renormalization
	exp_ftp1=strcat('i1/',mWM_ftp1);
	exp_ftp2=strcat('i1/',mWM_ftp2);

    clear matlabbatch
    
	% Run renormalization in ImCalc
	matlabbatch{1}.spm.util.imcalc.input = reg_ftp1_spm;
	matlabbatch{1}.spm.util.imcalc.output = renorm_ftp1;
	matlabbatch{1}.spm.util.imcalc.outdir = wd_spm;
	matlabbatch{1}.spm.util.imcalc.expression = exp_ftp1;
	matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
	matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
	matlabbatch{1}.spm.util.imcalc.options.mask = 0;
	matlabbatch{1}.spm.util.imcalc.options.interp = 1;
	matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
	matlabbatch{2}.spm.util.imcalc.input = reg_ftp2_spm;
	matlabbatch{2}.spm.util.imcalc.output = renorm_ftp2;
	matlabbatch{2}.spm.util.imcalc.outdir = wd_spm;
	matlabbatch{2}.spm.util.imcalc.expression = exp_ftp2;
	matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
	matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
	matlabbatch{2}.spm.util.imcalc.options.mask = 0;
	matlabbatch{2}.spm.util.imcalc.options.interp = 1;
	matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
	spm_jobman('run',matlabbatch);
    
end

if pet_mods==2||3

	% Reshape PET images for getting mean intensity
	PIB1=spm_vol(reg_pib1);
	PIB1_read=spm_read_vols(PIB1);
	rPIB1=reshape(PIB1_read,sz1*sz2*sz3,1);

	PIB2=spm_vol(reg_pib2);
	PIB2_read=spm_read_vols(PIB2);
	rPIB2=reshape(PIB2_read,sz1*sz2*sz3,1);

	% Get mean intensities
	ind_pib1=find(reroded_c2>0.5 & rPIB1>0);
	meanWM_pib1=mean(rPIB1(ind_pib1));
    mWM_pib1=num2str(meanWM_pib1);

	ind_pib2=find(reroded_c2>0.5 & rPIB2>0);
	meanWM_pib2=mean(rPIB2(ind_pib2));
    mWM_pib2=num2str(meanWM_pib2);

	% Define ImCalc expressions for renormalization
	exp_pib1=strcat('i1/',mWM_pib1);
	exp_pib2=strcat('i1/',mWM_pib2);
    
    clear matlabbatch

	% Run renormalization in ImCalc
	matlabbatch{1}.spm.util.imcalc.input = reg_pib1_spm;
	matlabbatch{1}.spm.util.imcalc.output = renorm_pib1;
	matlabbatch{1}.spm.util.imcalc.outdir = wd_spm;
	matlabbatch{1}.spm.util.imcalc.expression = exp_pib1;
	matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
	matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
	matlabbatch{1}.spm.util.imcalc.options.mask = 0;
	matlabbatch{1}.spm.util.imcalc.options.interp = 1;
	matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
	matlabbatch{2}.spm.util.imcalc.input = reg_pib2_spm;
	matlabbatch{2}.spm.util.imcalc.output = renorm_pib2;
	matlabbatch{2}.spm.util.imcalc.outdir = wd_spm;
	matlabbatch{2}.spm.util.imcalc.expression = exp_pib2;
	matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
	matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
	matlabbatch{2}.spm.util.imcalc.options.mask = 0;
	matlabbatch{2}.spm.util.imcalc.options.interp = 1;
	matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
	spm_jobman('run',matlabbatch);
    
end

clear matlabbatch

%Create brainmask
matlabbatch{1}.spm.util.imcalc.input = {avg_c1; avg_c2; avg_c3};
matlabbatch{1}.spm.util.imcalc.output = 'brainmask';
matlabbatch{1}.spm.util.imcalc.outdir = wd_spm;
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)>0.01';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);

% Define jacobians, brain mask, and output for masking of jacobians
[~,f1,e]=spm_fileparts(vols_mri(1,:));
[~,f2,~]=spm_fileparts(vols_mri(2,:));
jacobs=strcat(wd,'/jd_',f1,'_',f2,e);
mask=strcat(wd,'/brainmask.nii');
masked_output=strcat(wd,'/jd_',f1,'_',f2,'_masked');

clear matlabbatch

%Run ImCalc to mask jacobian images
matlabbatch{1}.spm.util.imcalc.input = {jacobs; mask};
matlabbatch{1}.spm.util.imcalc.output = masked_output;
matlabbatch{1}.spm.util.imcalc.outdir = wd_spm;
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
