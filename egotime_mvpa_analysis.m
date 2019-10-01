% Script for MVPA analysis based on the Princeton MVPA Toolbox
% Instead of using the 'VT_category-selective' mask of the
% tutorial, we will use a wholebrain mask that includes 43193 voxels.

% Adding the Path to Princeton-MVPA-toolbox
cd /home/vandamw/data2/MVPA_analyses/MATLAB/princeton-mvpa-toolbox-master
mvpa_add_paths;
addpath core

% initialize some variables
s= 1;
root = '/data2/user_data/vandamw/MVPA_analyses/MATLAB/ego_time/TS';
subjectvector = ['s05'];
subject = char(subjectvector(s,1:3));
direc=fullfile(root,subject);
direc2=fullfile(root,'MVPA_results');
output_dir = '/data2/user_data/vandamw/MVPA_analyses/MATLAB/ego_time/TS/MVPA_results';
cd (root)

cd (direc2)

% Here we create some csv files in which all of the results will be appended
fid = fopen('mvpa_accuracy_svm.csv','at')
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s','Run1','Run2','Run3','Run4','Run5','Run6');
fprintf(fid,'\n');
fclose(fid);

fid = fopen('mvpa_MVP_IVP_svm.csv','at')
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s','Run1','Run2','Run3','Run4','Run5','Run6');
fprintf(fid,'\n');
fclose(fid);

fid = fopen('mvpa_MW_IW_svm.csv','at')
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s','Run1','Run2','Run3','Run4','Run5','Run6');
fprintf(fid,'\n');
fclose(fid);

% At the end of the analysis we will append the data of each subject

cd (direc)

subj = init_subj('egotime','subject_s05')
subj = load_afni_mask(subj,'GreyMatterMask','GreyMatterMask+tlrc'); 
summarize(subj);

% We had 2 separate sessions with 4 runs each. We will train and test the classifier separately for
% each session.     
      
for i=1:6
    raw_filenames{i} = sprintf('mvpa_input_%s_run%i+tlrc',subject,i);
end

subj = load_afni_pattern(subj,'epi','GreyMatterMask',raw_filenames);

% At this point we have loaded all the epi data for the example subject
summarize(subj)
subj = init_object(subj,'regressors','conds');
load(sprintf('reg%s',subject));
subj = set_mat(subj,'regressors','conds',regs);
condnames = {'EC2','EC3','EM2','EM3','ST2','ST3','TC2','TC3','TM2','TM3'};
subj = set_objfield(subj,'regressors','conds','condnames',condnames);
condsobj = get_object(subj,'regressors','conds')

subj = init_object(subj,'selector','runs');
load('runs');
subj = set_mat(subj,'selector','runs',runs);
runs = get_object(subj,'selector','runs');
runs = get_mat(subj,'selector','runs');
summarize(subj);

% In this part of the script we will convolve the regressors with the HRF
% function. In order to do this we will need to call the
% convolve_regressors_afni.m file that is part of the toolbox. This program
% calls on the waver program of AFNI. However, it will need to know where
% this program is installed. It probably won't run if you just run the
% following command. I had to go into the convolve_regressors_afni.m file
% and insert the following statement 

% defaults.afni_location = '~/../Shared/abin'; # telling the program where
% AFNI is installed!


 %subj = convolve_regressors_afni(subj,'conds','runs', 'tr_size_in_seconds',1.85) % regressor
 
 %convolved with the HRF function. We, however want a binarized version!!
 
 %subj = shift_regressors(subj,'conds','runs',3)
 %summarize(subj,'objtype','regressors');

 subj = convolve_regressors_afni(subj,'conds','runs', ...
       'binarize_thresh',0.49, ...
       'tr_size_in_seconds',0.8, ...
       'do_plot', true);   % binarized version of the HRF convolved regressor
 summarize(subj,'objtype','regressors');
 conds_convt = get_mat(subj,'regressors','conds_convt')
 

%%
% First we will add the regressors for the 2 sentence presentations and 3
% sentence presentations, before adding them as a contrast to the 
% regressors we're interested in

% ----------------- A regressor consisting of the EgoMoving and the 
%                     EgoControl sentences

regs_EC_mat = get_mat(subj,'regressors','conds_convt');
regs_EC_mat([3:10],:) = [];
regs_EC = regs_EC_mat(1,:) + regs_EC_mat(2,:);
regs_EM_mat = get_mat(subj,'regressors','conds_convt');
regs_EM_mat([1:2 5:10],:) = [];
regs_EM = regs_EM_mat(1,:) + regs_EM_mat(2,:);
regs_EM_EC(1,:) = regs_EM;
regs_EM_EC(2,:) = regs_EC;
subj = init_object(subj,'regressors','regs_EM_EC');
subj = set_mat(subj,'regressors','regs_EM_EC', regs_EM_EC);

% ----------------- A regressor consisting of the TimeMoving and the 
%                     TimeControl sentences

regs_TC_mat = get_mat(subj,'regressors','conds_convt');
regs_TC_mat([1:6 9:10],:) = [];
regs_TC = regs_TC_mat(1,:) + regs_TC_mat(2,:);
regs_TM_mat = get_mat(subj,'regressors','conds_convt');
regs_TM_mat([1:8],:) = [];
regs_TM = regs_TM_mat(1,:) + regs_TM_mat(2,:);
regs_TM_TC(1,:) = regs_TM;
regs_TM_TC(2,:) = regs_TC;
subj = init_object(subj,'regressors','regs_TM_TC');
subj = set_mat(subj,'regressors','regs_TM_TC', regs_TM_TC);

% ------------------ A regressor consisting of the EgoMoving and the 
%                      TimeMoving sentences

regs_EM_TM(1,:) = regs_EM;
regs_EM_TM(2,:) = regs_TM;
subj = init_object(subj,'regressors','regs_EM_TM');
subj = set_mat(subj,'regressors','regs_EM_TM', regs_EM_TM);

% ------------------ A regressor consisting of the EgoMoving + TimeMoving
%                    versus the EgoControl + TimeControl sentences

regs_EMTM = regs_EM + regs_TM;
regs_ECTC = regs_EC + regs_TC;
regs_EMTM_ECTC(1,:) = regs_EMTM;
regs_EMTM_ECTC(2,:) = regs_ECTC;
subj = init_object(subj,'regressors','regs_EMTM_ECTC');
subj = set_mat(subj,'regressors','regs_EMTM_ECTC', regs_EMTM_ECTC);

% ------------------- A regressor consisting of the EgoControl and the 
%                       TimeControl sentences

regs_EC_TC(1,:) = regs_EC;
regs_EC_TC(2,:) = regs_TC;
subj = init_object(subj,'regressors','regs_EC_TC');
subj = set_mat(subj,'regressors','regs_EC_TC', regs_EC_TC);

% ------------------- A regressor consisting of the EgoMoving and the 
%                       Static Control sentences

regs_ST_mat = get_mat(subj,'regressors','conds_convt');
regs_ST_mat([1:4 7:10],:) = [];
regs_ST = regs_ST_mat(1,:) + regs_ST_mat(2,:);
regs_EM_ST(1,:) = regs_EM;
regs_EM_ST(2,:) = regs_ST;
subj = init_object(subj,'regressors','regs_EM_ST');
subj = set_mat(subj,'regressors','regs_EM_ST', regs_EM_ST);

% ------------------- A regressor consisting of the TimeMoving and the 
%                       Static Control sentences

regs_TM_ST(1,:) = regs_TM;
regs_TM_ST(2,:) = regs_ST;
subj = init_object(subj,'regressors','regs_TM_ST');
subj = set_mat(subj,'regressors','regs_TM_ST', regs_TM_ST);

% ------------------- A regressor consisting of EgoMoving + TimeMoving
%                       versus Static Control sentences

regs_EMTM_ST(1,:) = regs_EMTM;
regs_EMTM_ST(2,:) = regs_ST;
subj = init_object(subj,'regressors','regs_EMTM_ST');
subj = set_mat(subj,'regressors','regs_EMTM_ST', regs_EMTM_ST);
  
subj = create_norest_sel(subj,'conds_convt');
subj = create_norest_sel(subj,'regs_EM_EC');
subj = create_norest_sel(subj,'regs_TM_TC');
subj = create_norest_sel(subj,'regs_EM_TM');
subj = create_norest_sel(subj,'regs_EMTM_ECTC');
subj = create_norest_sel(subj,'regs_EC_TC');
subj = create_norest_sel(subj,'regs_EM_ST');
subj = create_norest_sel(subj,'regs_TM_ST');
subj = create_norest_sel(subj,'regs_EMTM_ST');

conds_convt_norest  = get_mat(subj,'selector','conds_convt_norest');
regs_EM_EC_norest   = get_mat(subj,'selector','regs_EM_EC_norest');
regs_TM_TC_norest   = get_mat(subj,'selector','regs_TM_TC_norest');
regs_EM_TM_norest   = get_mat(subj,'selector','regs_EM_TM_norest');
regs_EMTM_ECTC_norest = get_mat(subj,'selector','regs_EMTM_ECTC_norest');
regs_EC_TC_norest   = get_mat(subj,'selector','regs_EC_TC_norest');
regs_EM_ST_norest   = get_mat(subj,'selector','regs_EM_ST_norest');
regs_TM_ST_norest   = get_mat(subj,'selector','regs_TM_ST_norest');
regs_EMTM_ST_norest = get_mat(subj,'selector','regs_EMTM_ST_norest');

% here we will load the censor file that we created for each subject

subj = init_object(subj,'regressors','censor');
load(sprintf('censor%s.mat',subject));
subj = set_mat(subj,'regressors','censor',censor);

conds_convt_norest_outliersremoved = censor + conds_convt_norest;
regs_EM_EC_norest_outliersremoved = censor + regs_EM_EC_norest; 
regs_TM_TC_norest_outliersremoved = censor + regs_TM_TC_norest; 
regs_EM_TM_norest_outliersremoved = censor + regs_EM_TM_norest; 
regs_EMTM_ECTC_norest_outliersremoved = censor + regs_EMTM_ECTC_norest;
regs_EC_TC_norest_outliersremoved = censor + regs_EC_TC_norest;
regs_EM_ST_norest_outliersremoved = censor + regs_EM_ST_norest;
regs_TM_ST_norest_outliersremoved = censor + regs_TM_ST_norest;
regs_EMTM_ST_norest_outliersremoved = censor + regs_EMTM_ST_norest;

% creating a temporary regressor with all ones
regs_EM_EC_norest_outlier = ones(1,size(censor,2)); 
regs_TM_TC_norest_outlier = ones(1,size(censor,2));
regs_EM_TM_norest_outlier = ones(1,size(censor,2));
regs_EMTM_ECTC_norest_outlier = ones(1,size(censor,2));
regs_EC_TC_norest_outlier = ones(1,size(censor,2));
regs_EM_ST_norest_outlier = ones(1,size(censor,2)); 
regs_TM_ST_norest_outlier = ones(1,size(censor,2));
regs_EMTM_ST_norest_outlier = ones(1,size(censor,2)); 

% finding where either the outlier or norest regressor has a zero timepoint
regs_EM_EC_norest_outlier(find(regs_EM_EC_norest_outliersremoved~=2))= 0;
regs_TM_TC_norest_outlier(find(regs_TM_TC_norest_outliersremoved~=2))= 0;
regs_EM_TM_norest_outlier(find(regs_EM_TM_norest_outliersremoved~=2))= 0;
regs_EMTM_ECTC_norest_outlier(find(regs_EMTM_ECTC_norest_outliersremoved~=2))= 0;
regs_EC_TC_norest_outlier(find(regs_EC_TC_norest_outliersremoved~=2))= 0;
regs_EM_ST_norest_outlier(find(regs_EM_ST_norest_outliersremoved~=2))= 0; 
regs_TM_ST_norest_outlier(find(regs_TM_ST_norest_outliersremoved~=2))= 0;
regs_EMTM_ST_norest_outlier(find(regs_EMTM_ST_norest_outliersremoved~=2))= 0;

% now we will set the selector files for the newly created "sensor" files
% that are a combination of norest removed and outliers removed
subj = init_object(subj,'selector','regs_EM_EC_norest_outlier');
subj = set_mat(subj,'selector','regs_EM_EC_norest_outlier',regs_EM_EC_norest_outlier);
subj = init_object(subj,'selector','regs_TM_TC_norest_outlier');
subj = set_mat(subj,'selector','regs_TM_TC_norest_outlier',regs_TM_TC_norest_outlier);
subj = init_object(subj,'selector','regs_EM_TM_norest_outlier');
subj = set_mat(subj,'selector','regs_EM_TM_norest_outlier',regs_EM_TM_norest_outlier);
subj = init_object(subj,'selector','regs_EMTM_ECTC_norest_outlier');
subj = set_mat(subj,'selector','regs_EMTM_ECTC_norest_outlier',regs_EMTM_ECTC_norest_outlier);
subj = init_object(subj,'selector','regs_EC_TC_norest_outlier');
subj = set_mat(subj,'selector','regs_EC_TC_norest_outlier',regs_EC_TC_norest_outlier);
subj = init_object(subj,'selector','regs_EM_ST_norest_outlier');
subj = set_mat(subj,'selector','regs_EM_ST_norest_outlier',regs_EM_ST_norest_outlier);
subj = init_object(subj,'selector','regs_TM_ST_norest_outlier');
subj = set_mat(subj,'selector','regs_TM_ST_norest_outlier',regs_TM_ST_norest_outlier);
subj = init_object(subj,'selector','regs_EMTM_ST_norest_outlier');
subj = set_mat(subj,'selector','regs_EMTM_ST_norest_outlier',regs_EMTM_ST_norest_outlier);


%%
% zscoring the data by subtracting out the mean of each voxel's timecourse
% and scaling it so that the standard deviation of the timecourse is one.
subj = zscore_runs(subj,'epi','runs');
summarize(subj,'objtype','pattern');
before_zscore_mat = get_mat(subj,'pattern','epi');
after_zscore_mat = get_mat(subj,'pattern','epi_z');

%%
subj = init_object(subj,'regressors','conditions');
regs_conditions(1,:) = regs_EC;
regs_conditions(2,:) = regs_EM;
regs_conditions(3,:) = regs_ST;
regs_conditions(4,:) = regs_TC;
regs_conditions(5,:) = regs_TM;
subj = set_mat(subj,'regressors','conditions',regs_conditions);

 % In the following lines we are using the condition regressors and we are
 % creating blocklabels for each regressor so that we can add the
 % activations on a trial-by-trial basis
 
 %array of all ones

regressors_array_ones = ones(1,size(regs_conditions,2)); 
block_number = 0;
total_trials = length(regs_conditions);

regressors_array_ones(:,1:6) = 0;

for row = 1:size(regs_conditions,1)
for i = 1:(total_trials-7)
    if regs_conditions(row,i) == 0;
       regressors_array_ones(row,i+6) = 0;
    elseif regs_conditions(row,i) == 1 & regs_conditions(row,i-1) == 0;
       block_number = block_number + 1;
       regressors_array_ones(row,i+6) = block_number;
    elseif regs_conditions(row,i) == 1 & regs_conditions(row,i-1) == 1;
       regressors_array_ones(row,i+6) = block_number;            
    end
 
end

for i = total_trials-6
    if regs_conditions(row,i) == 0;
       regressors_array_ones(row,i+6) = 0;
    elseif regs_conditions(row,i) == 1 & regs_conditions(row,i-1) == 1;
       regressors_array_ones(row,i+6) = regressors_array_ones(row,i+5);
    end
end

end

all_regressors_array = regressors_array_ones(1,:) + regressors_array_ones(2,:) + regressors_array_ones(3,:) + regressors_array_ones(4,:) + regressors_array_ones(5,:);

%%
% Now we are adding the regressors_array_ones ARRAY to the subj structure
% as blocklabel variable, under the selector variables
subj = init_object(subj,'selector','blocklabels');
subj = set_mat(subj,'selector','blocklabels',all_regressors_array);

% here we add another variable which is the 'epi_z' pattern multiplied by
% the block variable.  
subj = average_object(subj,'regressors','conds_convt','blocklabels');
subj = average_object(subj,'regressors','regs_EM_EC','blocklabels');
subj = average_object(subj,'regressors','regs_TM_TC','blocklabels');
subj = average_object(subj,'regressors','regs_EM_TM','blocklabels');
subj = average_object(subj,'regressors','regs_EMTM_ECTC','blocklabels');
subj = average_object(subj,'regressors','regs_EC_TC','blocklabels');
subj = average_object(subj,'regressors','regs_EM_ST','blocklabels');
subj = average_object(subj,'regressors','regs_TM_ST','blocklabels');
subj = average_object(subj,'regressors','regs_EMTM_ST','blocklabels');
subj = average_object(subj,'pattern','epi_z','blocklabels');

% Here we will create a regressor to do a crossmodal classification
% basically we just add the MVP,MW conditions and the IW and IVP conditions

regs_MVP_IVP_avg = get_mat(subj,'regressors','regs_MVP_IVP_avg');
regs_MW_IW_avg = get_mat(subj,'regressors','regs_MW_IW_avg');
regs_crossmodal_avg = regs_MVP_IVP_avg + regs_MW_IW_avg;
subj = init_object(subj,'regressors','regs_crossmodal_avg');
subj = set_mat(subj, 'regressors','regs_crossmodal_avg', regs_crossmodal_avg);

creating the group of cross-validation selectors, train or test on all
pictures or all words

regs_MVP_IVP_avg_ones = regs_MVP_IVP_avg(1,:) + regs_MVP_IVP_avg(2,:);
regs_MVP_IVP_avg_threes = (regs_MVP_IVP_avg(1,:) + regs_MVP_IVP_avg(2,:))*3;
regs_MW_IW_avg_ones = regs_MW_IW_avg(1,:) + regs_MW_IW_avg(2,:);
regs_MW_IW_avg_threes = (regs_MW_IW_avg(1,:) + regs_MW_IW_avg(2,:))*3;

%%
% creating the group of cross-validation selectors, there will be one for
% each iteration of the cross-validation process. 

runs = get_mat(subj,'selector','runs');
nRuns = max(runs);
nTimepoints = length(runs);

runs_xval_sl = ones(nRuns,nTimepoints);

for r =1:nRuns
    cur_final_test_run = find(runs==r);
    runs_xval_sl(r, cur_final_test_run) = 1; %normally here you put your 2s, in our case we don't do any feature selection!
end

imagesc(runs_xval_sl)
set(gca,'Clim',[1 3])
colorbar

% At this point we will add the searchlight-generalization timepoints (3s
% in the matrix

for r=1:nRuns
    cur_searchlight_gen_run = find(runs== nRuns-r+1);
    runs_xval_sl(r, cur_searchlight_gen_run) = 3;
end

imagesc(runs_xval_sl)
colorbar

%%
% In this part we will multiply the runs_xval_sl selector with the no_rest_outlier censor files. This ensures
% that we only train and test on datapoints that have one of the conditions of interest.

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EM_EC_norest_outlier;
    runs_norest_xval_EM_EC(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_TM_TC_norest_outlier;
    runs_norest_xval_TM_TC(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EM_TM_norest_outlier;
    runs_norest_xval_EM_TM(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EMTM_ECTC_norest_outlier;
    runs_norest_xval_EMTM_ECTC(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EC_TC_norest_outlier;
    runs_norest_xval_EC_TC(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EM_ST_norest_outlier;
    runs_norest_xval_EM_ST(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_TM_ST_norest_outlier;
    runs_norest_xval_TM_ST(i,:) = a(i,:).*b(i,:);
end

for i = 1:nRuns
    a(i,:) = runs_xval_sl(i,:);
    b(i,:) = regs_EMTM_ST_norest_outlier;
    runs_norest_xval_EMTM_ST(i,:) = a(i,:).*b(i,:);
end


%%
% here the runs_norest_xval files are added to the subject structure

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_EC_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EM_EC(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_TM_TC_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_TM_TC(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_TM_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EM_TM(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EMTM_ECTC_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EMTM_ECTC(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EC_TC_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EC_TC(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_ST_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EM_ST(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_TM_ST_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_TM_ST(r,:));
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EMTM_ST_%i',r);
    subj = initset_object(subj,'selector', cur_name, ...
                          runs_norest_xval_EMTM_ST(r,:));
end

%%
% At this point we will multiply the selectors with the blocklabel variable
% and we remove the original runs_norest_xval files

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_EC_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_TM_TC_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_TM_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EMTM_ECTC_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EC_TC_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EM_ST_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_TM_ST_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

for r=1:nRuns
    cur_name = sprintf('runs_norest_xval_EMTM_ST_%i',r);
    subj = average_object(subj,'selector', cur_name,'blocklabels');
    subj = remove_object(subj,'selector',cur_name);
end

%%
% Here we will get the newly created runs_norest_xval_avg files and will
% add them to the subj structure under a grouping variable

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EM_EC_%i_avg',r);
     runs_norest_xval_EM_EC_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EM_EC_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EM_EC_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_TM_TC_%i_avg',r);
     runs_norest_xval_TM_TC_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_TM_TC_avg(r,:), ...
                           'group_name', 'runs_norest_xval_TM_TC_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EM_TM_%i_avg',r);
     runs_norest_xval_EM_TM_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EM_TM_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EM_TM_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EMTM_ECTC_%i_avg',r);
     runs_norest_xval_EMTM_ECTC_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EMTM_ECTC_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EMTM_ECTC_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EC_TC_%i_avg',r);
     runs_norest_xval_EC_TC_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EC_TC_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EC_TC_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EM_ST_%i_avg',r);
     runs_norest_xval_EM_ST_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EM_ST_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EM_ST_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_TM_ST_%i_avg',r);
     runs_norest_xval_TM_ST_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_TM_ST_avg(r,:), ...
                           'group_name', 'runs_norest_xval_TM_ST_avg');

end

for r=1:nRuns
     cur_name = sprintf('runs_norest_xval_EMTM_ST_%i_avg',r);
     runs_norest_xval_EMTM_ST_avg(r,:) = get_mat(subj,'selector', cur_name);
     subj = remove_object(subj,'selector',cur_name);
     subj = initset_object1(subj,'selector', cur_name, ...
                           runs_norest_xval_EMTM_ST_avg(r,:), ...
                           'group_name', 'runs_norest_xval_EMTM_ST_avg');

end

%%

% In this part we will create the adjacency matrix
subj.adj_sphere = create_adj_list(subj,'GreyMatterMask','radius',3); %the default radius is 2.0 we can change this by setting RAADIUS
%subj.adj_sphere = create_adj_list(subj,'wholebrain','radius',1); %the default radius is 2.0 we can change this by setting RAADIUS

%%

class_args.train_funct_name = 'train_linearsvmlog';
class_args.test_funct_name = 'test_linearsvmlog';
class_args.nHidden = 0;

scratch.class_args = class_args;
scratch.perfmet_funct = 'perfmet_maxclass';
scratch.perfmet_args = struct([]);

statmap_srch_arg.adj_list = subj.adj_sphere;
statmap_srch_arg.obj_funct = 'statmap_classify';
statmap_srch_arg.scratch = scratch;

%%

subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EM_EC_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EM_EC_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EM_EC_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_EC_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_EC_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_TM_TC_avg', ... % binary regs (for GNB)
       'runs_norest_xval_TM_TC_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_TM_TC_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_TC_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_TC_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EM_TM_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EM_TM_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EM_TM_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_TM_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_TM_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EMTM_ECTC_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EMTM_ECTC_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EMTM_ECTC_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ECTC_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ECTC_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
    
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EC_TC_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EC_TC_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EC_TC_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EC_TC_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EC_TC_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EM_ST_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EM_ST_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EM_ST_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EM_ST_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EM_ST_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_TM_ST_avg', ... % binary regs (for GNB)
       'runs_norest_xval_TM_ST_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_TM_ST_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_TM_ST_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/TM_ST_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
subj = feature_select( ...
       subj, ...
       'epi_z_avg', ... % data
       'regs_EMTM_ST_avg', ... % binary regs (for GNB)
       'runs_norest_xval_EMTM_ST_avg', ... % selector      %replaced the 'runs_xval_sl' selector from the example!!
       'statmap_funct','statmap_searchlight', ...% function
       'statmap_arg',statmap_srch_arg, ...
       'new_map_patname', 'epi_z_srch_EMTM_ST_avg_no_sel', ...
       'thresh',[]);
   
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_1','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_1_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_2','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_2_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_3','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_3_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_4','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_4_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_5','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_5_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_6','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_6_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_7','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_7_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
write_to_afni(subj,'pattern','epi_z_srch_EMTM_ST_avg_no_sel_8','GreyMatterMask+tlrc','output_filename',sprintf('%s/EMTM_ST_allvoxels_8_searchlight_results_svm_avg_no_sel_%s',output_dir,subject))
  
clear all


