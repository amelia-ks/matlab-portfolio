%%----------------------------------------------------------------------------
% Pipeline to fully process Florbetapir (FBP) scans according to ADNI process.
% Output is FBP SUVR image, extraction and reference region masks, and SUVR
% quantification (and most files created along the way).
%
% Adapted from a defunct Python script by Amelia, May 2020
%
% To do before running:
%   1. Run all subjects through Freesurfer, make sure results are available
%      in directory named "fs_dir" below
%   2. Transfer dicoms to raw frames directory and run the AV45_Conversion
%      script to create 4 nifti frames
%   3. Edit china_basin_pet_tracking.xlsx spreadsheet to include scan/MRI
%      info for each subject
%%----------------------------------------------------------------------------

% Define relevant directories
proj_dir='/home/jagust/UCSF/china_basin_pet/AV45/';
scans_dir=strcat(proj_dir,'scans/');
raw_dir=strcat(scans_dir,'raw_frames/');
fs_dir='/home/jagust/UCSF/freesurfer_5_1-v2/';

% Change subjects to be run here!
subs=char('id1', 'id2');

T=readtable('/home/jagust/UCSF/china_basin_pet/AV45/info/china_basin_pet_tracking.xlsx'); %Load tracking spreadsheet
ppidnT=strcat('p',T.pidn); %Add "p" to all PIDNs

spm('defaults','PET');

% Find number of subjects to loop through
numbersubs=size(subs,1);

% Preallocate fields that will be filled in the loop
av45date{numbersubs,1}=[];
mridate{numbersubs,1}=[];
cingulate_mean=zeros(size(subs,1),1);
frontal_mean=zeros(size(subs,1),1);
parietal_mean=zeros(size(subs,1),1);
temporal_mean=zeros(size(subs,1),1);
suvr_quant=zeros(size(subs,1),1);

%Loop through subjects for processing
for i=1:size(subs,1)

    sub=subs(i,:); %Define individual subject
    subind=find(strcmp(cellstr(ppidnT),cellstr(sub))); %Find row associated with subject in spreadsheet

    % Find subject-specific directories
    sub_rawdir=strcat(raw_dir,sub); %for raw frames
    sub_fsdir=strcat(fs_dir,sub); %for mri files

    % Define scan/mri dates and create strings
    % Need to add 693960 to each date number to be accurate...
    d1temp=datenum(T.av45_date(subind))+693960;
    av45dateS=datestr(d1temp,'mm-dd-yyyy');
    av45date(i,1)=cellstr(av45dateS); %to use in final table
    d2temp=datenum(T.mri_date(subind))+693960;
    mridateS=datestr(d2temp,'mm-dd-yyyy');
    mridate(i,1)=cellstr(mridateS); %to use in final table

    % Create subject-specific scans directory
    sub_dir=strcat(scans_dir,sub);
    mkdir (sub_dir)

    % Create directories for processing output
    sub_av45dir=strcat(sub_dir,'/av45_',av45dateS);
    sub_mridir=strcat(sub_dir,'/mri_',mridateS);
    mkdir (sub_av45dir)
    mkdir (sub_mridir)
    sub_realigndir=strcat(sub_av45dir,'/realign');
    mkdir (sub_realigndir)

    % Find and count frames, copy over to realign directory
    frame_path=strcat(sub_rawdir,'/NEWoutput*nii');
    framefiles=dir(frame_path);
    num_frames=size(framefiles,1); %counts number of frames
    copyfile (frame_path,sub_realigndir)

    % Delete subject-specific raw frame directory
    allraw=strcat(sub_rawdir,'/*');
    delete(allraw)
    rmdir(sub_rawdir)

    % Preallocate frame cell arrays
    frame{num_frames,1}=[];
    frame_spm{num_frames,1}=[];

    % Rename frames, assuming there are 4
    for n=1:num_frames

        frame_old=strcat(sub_realigndir,'/NEWoutput',num2str(n),'.nii'); %path now
        frame(n)=cellstr(strcat(sub_realigndir,'/',sub,'_av45_fr',num2str(n),'.nii')); %desired path
        frame_spm(n)=cellstr(strcat(char(frame(n)),',1')); %Prepping frames for SPM language
        movefile (frame_old,char(frame(n)))

    end

    % Find nu.mgz and aparc+aseg.mgz files in FS directory and copy over
    numgz_path=strcat(sub_fsdir,'/mri/nu.mgz');
    aparcmgz_path=strcat(sub_fsdir,'/mri/aparc+aseg.mgz');
    copyfile (numgz_path,sub_mridir)
    copyfile (aparcmgz_path,sub_mridir)

    % Define mgz/nii aparc and nu files for easy access
    numgz=strcat(sub_mridir,'/nu.mgz');
    aparcmgz=strcat(sub_mridir,'/aparc+aseg.mgz');
    aparc=strcat(sub_mridir,'/aparc+aseg.nii');
    nu=strcat(sub_mridir,'/nu.nii');

    % Create strings to convert mgz to nii, use Linux system to do it
    convertnu=char(strcat('mri_convert',{' '},numgz,{' '},nu));
    convertaparc=char(strcat('mri_convert',{' '},aparcmgz,{' '},aparc));
    system(convertnu); %deploys Linux
    system(convertaparc); %deploys Linux

    % Remove mgz files
    delete (numgz)
    delete (aparcmgz)

    clear smoothbatch

    % Smooth frames with a [5.8 x 5.8 x 5.84] kernel
    smoothbatch{1}.spm.spatial.smooth.data = frame_spm;
    smoothbatch{1}.spm.spatial.smooth.fwhm = [5.8 5.8 5.84];
    smoothbatch{1}.spm.spatial.smooth.dtype = 0;
    smoothbatch{1}.spm.spatial.smooth.im = 0;
    smoothbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',smoothbatch);

    % Preallocate smoothed frame array
    sframe{num_frames,1}=[];

    % Define smoothed frames for realign batch, delete unsmoothed frames
    for m=1:num_frames

        sframe(m)=cellstr(strcat(sub_realigndir,'/s',sub,'_av45_fr',num2str(m),'.nii,1'));
        frameUS=char(frame(m));
        delete (frameUS)

    end

    clear realignbatch

    % Realign smoothed frames and create mean image
    realignbatch{1}.spm.spatial.realign.estwrite.data = {sframe};
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realignbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    realignbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realignbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realignbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realignbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realignbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',realignbatch);

    % Define path to mean image
    mean_path=strcat(sub_realigndir,'/means',sub,'_av45_fr1.nii');
    mean_spm=cellstr(strcat(sub_mridir,'/means',sub,'_av45_fr1.nii,1')); %for spm

    % Copy mean image to mri path
    copyfile (mean_path,sub_mridir)

    % Define aparc and nu nii files for easy access
    nu_spm=strcat(nu,',1'); %for spm
    aparc_spm=strcat(aparc,',1'); %for spm

    clear coregbatch

    % Coregister mean image to nu image
    coregbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(nu_spm);
    coregbatch{1}.spm.spatial.coreg.estwrite.source = mean_spm;
    coregbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    coregbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    coregbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    coregbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    coregbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    coregbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4; %4th deg
    coregbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    coregbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    coregbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run',coregbatch);

    clear resliceaparcbatch

    % Reslice aparc to nu to ensure voxel sizes/coordinate systems are
    % consistent
    resliceaparcbatch{1}.spm.spatial.coreg.write.ref = cellstr(nu_spm);
    resliceaparcbatch{1}.spm.spatial.coreg.write.source = cellstr(aparc_spm);
    resliceaparcbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; %NN
    resliceaparcbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    resliceaparcbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    resliceaparcbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',resliceaparcbatch);

    % Define path to coregistered mean and resliced aparc
    rmean=strcat(sub_mridir,'/rmeans',sub,'_av45_fr1.nii,1');
    raparc=strcat(sub_mridir,'/raparc+aseg.nii');
    raparc_spm=strcat(raparc,',1');

    % Delete non-coregistered mean from mri folder
    mean_mripath=strcat(sub_mridir,'/means',sub,'_av45_fr1.nii');
    delete(mean_mripath)

    % Remove non-resliced aparc file
    delete (aparc)

    % Combine aparc and rmean to be used in creating cerebellum mask
    tempvol=[cellstr(raparc_spm);cellstr(rmean)];

    clear cermaskbatch

    % Create mask of whole cerebellum as reference region
    cermaskbatch{1}.spm.util.imcalc.input = tempvol;
    cermaskbatch{1}.spm.util.imcalc.output = 'whole_cereb_mask';
    cermaskbatch{1}.spm.util.imcalc.outdir = cellstr(sub_mridir);
    cermaskbatch{1}.spm.util.imcalc.expression = '(((i1==8)+(i1==47)+(i1==7)+(i1==46))>0).*(i2>0)';
    cermaskbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    cermaskbatch{1}.spm.util.imcalc.options.dmtx = 0;
    cermaskbatch{1}.spm.util.imcalc.options.mask = 0;
    cermaskbatch{1}.spm.util.imcalc.options.interp = 0; %NN
    cermaskbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',cermaskbatch);

    % Define path to cerebellar mask
    cereb_mask=strcat(sub_mridir,'/whole_cereb_mask.nii');
    cereb_mask_spm=strcat(cereb_mask,',1'); %for spm

    clear reslicemaskbatch

    % Reslice cerebellar mask to ensure it shares voxel size/coordinate
    % system with rmean image
    reslicemaskbatch{1}.spm.spatial.coreg.write.ref = cellstr(rmean);
    reslicemaskbatch{1}.spm.spatial.coreg.write.source = cellstr(cereb_mask_spm);
    reslicemaskbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; %NN
    reslicemaskbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    reslicemaskbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    reslicemaskbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',reslicemaskbatch);

    % Define resliced cerebellar mask
    rcer_mask=strcat(sub_mridir,'/rwhole_cereb_mask.nii');
    rcermask_spm=strcat(rcer_mask,',1');

    % Remove non-resliced cerebellar mask
    delete (cereb_mask)

    % Extract info from cerebellar mask to extract values
    Cere1=spm_vol(rcer_mask);
	Cere2=spm_read_vols(Cere1);

    % Read and resize mean image to extract values
    Mean1=spm_vol(rmean); %loads header info
	Mean2=spm_read_vols(Mean1); %loads images of interest into FDG2
	Mean_mask=Mean2.*Cere2; %mask mean PET with cerebellar mask
    Mean_mask=nonzeros(Mean_mask); %remove any 0 voxels
    Mean_mask=Mean_mask(~isnan(Mean_mask)); %remove any 0 voxels
    masksz=size(Mean_mask,1); %determine cluster size of reference ROI

    % Get mean within reference region
    ref_mean=mean(Mean_mask);

    % Display reference region info for user
    fprintf('\n\n**\nMean in whole cerebellum for %s is %d, extracted from %d voxels. \nI will now create the SUVR using this value!\n**\n\n',sub,ref_mean,masksz)

    clear suvrbatch

    % Divide coregistered mean by reference region mean
    exp=strcat('i1/',num2str(ref_mean)); %define imcalc expression
    output=strcat(sub,'_av45_rsuvr_wholecer'); %name suvr output

    % Create SUVR image using reference region mean
    suvrbatch{1}.spm.util.imcalc.input = cellstr(rmean);
    suvrbatch{1}.spm.util.imcalc.output = output;
    suvrbatch{1}.spm.util.imcalc.outdir = cellstr(sub_mridir); %will output SUVR to mri directory
    suvrbatch{1}.spm.util.imcalc.expression = exp;
    suvrbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    suvrbatch{1}.spm.util.imcalc.options.dmtx = 0;
    suvrbatch{1}.spm.util.imcalc.options.mask = 0;
    suvrbatch{1}.spm.util.imcalc.options.interp = 1; %trilinear
    suvrbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',suvrbatch);

    % Define SUVR path
    suvr=strcat(sub_mridir,'/',output,'.nii');
    suvr_spm=strcat(sub_mridir,'/',output,'.nii,1'); %for spm
    tempvol2=[cellstr(raparc_spm);cellstr(suvr_spm)]; %combine for batch

    % Invert suvr for clinical reads using FSL in Linux
    invert=char(strcat('fslmaths',{' '},suvr,{' '},'-div',{' '},'-1',{' '},sub_mridir,'/',output,'_inverted'));
    system(invert);
    gunzip=char(strcat('gunzip',{' '},sub_mridir,'/',output,'_inverted.nii.gz'));
    system(gunzip);

    clear extmaskbatch

    % Create extraction ROI mask for easy visualization, includes all lobes
    extmaskbatch{1}.spm.util.imcalc.input = tempvol2;
    extmaskbatch{1}.spm.util.imcalc.output = 'extraction_ROI';
    extmaskbatch{1}.spm.util.imcalc.outdir = cellstr(sub_mridir);
    extmaskbatch{1}.spm.util.imcalc.expression = '(((i1==1015)+(i1==1030)+(i1==2015)+(i1==2030)+(i1==1008)+(i1==1025)+(i1==1029)+(i1==1031)+(i1==2008)+(i1==2025)+(i1==2029)+(i1==2031)+(i1==1003)+(i1==1032)+(i1==1012)+(i1==1014)+(i1==1018)+(i1==1019)+(i1==1020)+(i1==1027)+(i1==1028)+(i1==2003)+(i1==2032)+(i1==2012)+(i1==2014)+(i1==2018)+(i1==2019)+(i1==2020)+(i1==2027)+(i1==2028)+(i1==1002)+(i1==1010)+(i1==1023)+(i1==1026)+(i1==2002)+(i1==2010)+(i1==2023)+(i1==2026))>0).*(i2>0)';
    extmaskbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    extmaskbatch{1}.spm.util.imcalc.options.dmtx = 0;
    extmaskbatch{1}.spm.util.imcalc.options.mask = 0;
    extmaskbatch{1}.spm.util.imcalc.options.interp = 1;
    extmaskbatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',extmaskbatch);

    % Define extraction ROI for use
    extmaskroi=strcat(sub_mridir,'/extraction_ROI.nii');
    extmaskroi_spm=strcat(sub_mridir,'/extraction_ROI.nii,1');

    % Read and resize suvr image to extract lobar means
    SUVR1=spm_vol(suvr);
    SUVR2=spm_read_vols(SUVR1);
    [sz1,sz2,sz3]=size(SUVR2);
    SUVR3=reshape(SUVR2,sz1*sz2*sz3,1);

    % Read and resize aparc image to extract indices
    Aparc1=spm_vol(raparc);
    Aparc2=spm_read_vols(Aparc1);
    Aparc3=reshape(Aparc2,sz1*sz2*sz3,1);

    % Define lobar ROIs
    cing_index=find(Aparc3==1002 | Aparc3==1010 | Aparc3==1023 | Aparc3==1026 | Aparc3==2002 | Aparc3==2010 | Aparc3==2023 | Aparc3==2026 & SUVR3>0);
    front_index=find(Aparc3==1003 | Aparc3==1032 | Aparc3==1012 | Aparc3==1014 | Aparc3==1018 | Aparc3==1019 | Aparc3==1020 | Aparc3==1027 | Aparc3==1028 | Aparc3==2003 | Aparc3==2032 | Aparc3==2012 | Aparc3==2014 | Aparc3==2018 | Aparc3==2019 | Aparc3==2020 | Aparc3==2027 | Aparc3==2028 & SUVR3>0);
    parietal_index=find(Aparc3==1008 | Aparc3==1025 | Aparc3==1029 | Aparc3==1031 | Aparc3==2008 | Aparc3==2025 | Aparc3==2029 | Aparc3==2031 & SUVR3>0);
    temporal_index=find(Aparc3==1015 | Aparc3==1030 | Aparc3==2015 | Aparc3==2030 & SUVR3>0);

    % Get means in each lobe
    cingulate_mean(i,1)=mean(SUVR3(cing_index));
    frontal_mean(i,1)=mean(SUVR3(front_index));
    parietal_mean(i,1)=mean(SUVR3(parietal_index));
    temporal_mean(i,1)=mean(SUVR3(temporal_index));

    % AV45  quant=sum(4 lobe regions)/4; i.e. non-weighted mean
    suvr_quant(i,1)=(cingulate_mean(i,1)+frontal_mean(i,1)+parietal_mean(i,1)+temporal_mean(i,1))/4;

    % Combine both masks to create combined ROI for visualization
    masks=[cellstr(rcermask_spm);cellstr(extmaskroi_spm)];

    clear masksroibatch

    % Combine reference region and extraction ROI in same image for
    % visualization
    masksroibatch{1}.spm.util.imcalc.input = masks;
    masksroibatch{1}.spm.util.imcalc.output = 'reference_extraction_ROIs';
    masksroibatch{1}.spm.util.imcalc.outdir = cellstr(sub_mridir);
    masksroibatch{1}.spm.util.imcalc.expression = '(i1*1)+(i2*2)';
    masksroibatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    masksroibatch{1}.spm.util.imcalc.options.dmtx = 0;
    masksroibatch{1}.spm.util.imcalc.options.mask = 0;
    masksroibatch{1}.spm.util.imcalc.options.interp = 1;
    masksroibatch{1}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',masksroibatch);

    % Define image
    bothmasks=strcat(sub_mridir,'/reference_extraction_ROIs.nii');

    % Remove image with just extraction ROI, unnecessary
    delete (extmaskroi)

    % Display quant for user
    fprintf('\n\n**\nCompleted processing for %s.\nSUVR quantification is %s.\n**\n\n',sub,num2str(suvr_quant(i,1)));

end

% Prep values to display in results file
C=[cellstr(subs),av45date,mridate,num2cell(cingulate_mean),num2cell(frontal_mean),num2cell(parietal_mean),num2cell(temporal_mean),num2cell(suvr_quant)];
FinalT=cell2table(C);
resultsfile='/home/jagust/UCSF/china_basin_pet/AV45/info/av45_suvr_values.csv'; %results found here
FinalT.Properties.VariableNames={'Subject' 'AV45_Date' 'MRI_Date' 'Cingulate_Mean' 'Frontal_Mean' 'Lateral_Parietal_Mean' 'Lateral_Temporal_Mean' 'AV45_Index'};
writetable(FinalT,resultsfile)

clear;
