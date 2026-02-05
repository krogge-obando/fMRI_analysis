%This code will derive the whole-brain seedbased correaltion 
% of the insula for all the parklab data.


%derive masks, and rois
maskfile = '/data/neurogroup/parklab/MNI152_T1_2mm_brain_Mask.nii';
mm = niftiread(maskfile);
bm = find(mm==1);

insula_roi=niftiread("/data/neurogroup/parklab/ROIs/INSULA_ROI_AAL3.nii");

insula_vox=find(insula_roi>0);

[~,~,ir_insula] = intersect(insula_vox,bm);

% Now we will prep the data to loop through the subjects

proc_dir = '/data/neurogroup/parklab_tmp/';

sub_list= dir([proc_dir,'sub*']);
%%
% write the correlation maps (all together)
corrmaps_dir = '/data/neurogroup/parklab/PPS_study/ROI_maps/insula_2326/';
mkdir(corrmaps_dir)  
% write Ymat files
Ymat_dir = '/data/neurogroup/parklab/PPS_study/Ymats/insula_2326';
mkdir(Ymat_dir);

for qq=5:32
       
    sub_id = sub_list(qq).name
    
    % single & multi-echo scans
    se_scans = dir([proc_dir,sub_id,'/proc_scan*']);
    me_scans = dir([proc_dir,sub_id,'/meica_proc_scan*']);  
    scans = cat(1,se_scans,me_scans);
    
    for s=1:length(scans)
        scan_id = scans(s).name ;
        if strcmp(scan_id(1:3),'mei');
            in_dir = [proc_dir,sub_id,'/',scan_id,'/meica_out/'];
        else
            in_dir = [proc_dir,sub_id,'/',scan_id,'/'];
        end
        scan_no = scan_id(end-5:end);
        scan_prefix = [sub_id,'-',scan_no];

        %check if fMRI exists

        fmri_file = [in_dir,scan_prefix,'_EPI2MNI_sm_nr.nii.gz'];
        
        
        if ~isfile(fmri_file)

            fprintf('File "%s" not found.\n',fmri_file)

            continue;

        end



        fprintf('File "%s" found. Proceeding with code.', fmri_file);

       
        nn = niftiread(fmri_file);
        dims=size(nn);
        V = reshape(nn, prod(dims(1:3)),[]);
        Y = V(bm,:)';
        
        % extract seed time courses
       
        tc_insula = mean(Y(:,ir_insula),2);
        
        % correlation map
         
        map_insula_r = corr(tc_insula,Y);

        % fisher z
        
        map_insula_z = atanh(map_insula_r);
        

        % vols
        
        vol_insula_z = zeros(91,109,91);
        vol_insula_z(bm) = map_insula_z;
        
        
        masknii = niftiread(maskfile);
        nii_template = masknii;
        template_info = niftiinfo(maskfile);
        template_info.Datatype = 'double';
        
        %derive as a nifti file
       
        niftiwrite(vol_insula_z,[corrmaps_dir,scan_prefix,'_INSULA_z.nii'],template_info);


        % also save Ymat in a folder (all together)
    outname = [Ymat_dir,'/',scan_prefix,'_Ymni_sm_nr.mat'];
    sourcefile = fmri_file;
    display(['saving as: ',outname]);
    save(outname,'Y','sourcefile','maskfile');
     end 
    clear Y tc_insula 
    end






