function trans_full = trans_prob(rangeID,FreqRange)
numS = length(FreqRange); 

Training_data_dir = '..\..\datasets\Training_data\';
Extra_TrainingData_dir =  '..\..\datasets\Extra_TrainingData\';
TestData_dir = '..\..\datasets\TestData\';
TrueBPM_dir = '..\..\datasets\TrueBPM\';

IDData = {'DATA_01_TYPE01','DATA_02_TYPE02','DATA_03_TYPE02','DATA_04_TYPE02',...
    'DATA_05_TYPE02','DATA_06_TYPE02','DATA_07_TYPE02','DATA_08_TYPE02','DATA_09_TYPE02',...
    'DATA_10_TYPE02','DATA_11_TYPE02','DATA_12_TYPE02','DATA_S04_T01',...
    'TEST_S01_T01', 'TEST_S02_T01', 'TEST_S02_T02', 'TEST_S03_T02', ...
    'TEST_S04_T02', 'TEST_S05_T02', 'TEST_S06_T01', 'TEST_S06_T02',...
    'TEST_S07_T02', 'TEST_S08_T01'};


trans_full = zeros(numS,numS);
for idnb=rangeID
    if idnb<= 12
        load([Training_data_dir IDData{idnb} '_BPMtrace']);           % load groundtruth
    elseif idnb == 13
        load([Extra_TrainingData_dir 'BPM_S04_T01.mat']);           % load groundtruth
    else
        load([TrueBPM_dir 'True' IDData{idnb}(5:end)]);           % load groundtruth
    end
    %fullBPM0 = [ fullBPM0 BPM0' ];
    
    %# discretization/quantization into 8 levels
    x = BPM0;
    edges = linspace(60,180,numS+1);
    [counts,bins] = histc(x, edges);
    
    bins(bins==0) = 1;
    counts(1) = sum(bins==1);
    %# fix last level of histc output
    last = numel(counts);
    bins(bins==last) = last - 1;
    counts(last-1) = counts(last-1) + counts(last);
    counts(last) = [];
    
    %# show histogram
    %bar(edges(1:end-1), counts, 'histc')
    
    %# transition matrix
    trans = full(sparse(bins(1:end-1), bins(2:end), 1));
    dim = setdiff(1:numS,1:size(trans,2));
    trans(dim,dim)=0;
    trans_full = trans_full + trans;
end
trans_full = bsxfun(@rdivide, trans_full, sum(trans_full,2));
trans_full(isnan(trans_full)) = eps;
