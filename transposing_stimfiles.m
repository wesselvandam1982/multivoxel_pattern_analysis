%MATLAB script for transposing stimfiles
subject_list = {'203','204','205','206','207','208','209','210','211','212','213','214','216','217','218','219','220'}
k = 203;

for s = 1:length(subject_list)
    
pwd
cd (subject_list{s})

list = dir('*.1D')
for i = 1:length(list)
    file{i} = getfield(list,{i},'name')
    file_read{i} = textread(file{i});
    file_trans{i} = transpose(file_read{i});
    k = str2num(subject_list{s})
end
    
    regs = [file_trans{1};file_trans{2};file_trans{3};file_trans{4};file_trans{5};file_trans{6};file_trans{7};file_trans{8};file_trans{9}]
    fname = sprintf('regs%d.mat',(k));
    save (fname,'regs')
    clear file
    clear file_read
    clear file_trans

%clear all

cd ..

end