%MATLAB script for creating outlier censor mat files

subject_list = {'210','211','212','213','214','216','217','218','219','220'}
%k = 203;

for s = 1:length(subject_list)
    
pwd
cd (subject_list{s})
k = str2num(subject_list{s})
folder = sprintf('sCSFe_tshift_outlier_procedure_MVPA_input%d.results',(k));
cd (folder)
pwd

list = dir('*combined_2.1D')      %here it graps the censor file that is outputted by AFNI that will tell which datapoints are censored out. 
for i = 1:length(list)
    file{i} = getfield(list,{i},'name')
    file_read{i} = textread(file{i});
    file_trans{i} = transpose(file_read{i});
    k = str2num(subject_list{s})
end
    
    censor_2512points = [file_trans{1}];
    censor = censor_2512points([1],[1:310 315:624 629:938 943:1252 1257:1566 1571:1880 1885:2194 2199:2508]); % there is 4 extra volumes at the end of each run we will have to get rid off
    fname = sprintf('censor%d.mat',(k));
    save (fname,'censor')
    movefile censor*.mat ..
    clear file
    clear file_read
    clear file_trans

%clear all

cd ../..

end
