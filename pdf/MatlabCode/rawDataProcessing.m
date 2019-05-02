clear

%Convert and align raw data input
%See documentation for suggested data format

%30d implied volatility
vol30dATM = readtable('SPX Vol 30d.csv');
vol30dATM = vol30dATM(vol30dATM.delta==-50,:);
vol30dATM = vol30dATM(vol30dATM.days==30,:);

%option data
OptionData_RAW = readtable('SPX Option Data 1996-2017 Full Maturity.csv');
%underlying price
IndexData = readtable('SPX Index 1996-2017.csv');

OptionData_Clean = zeros(12);
IndexData_Clean = zeros(2);

lastDate=1;
IndexData = [datenum(int2str(IndexData{:,2}),'yyyymmdd') IndexData{:,7}];

parfor i = 1:height(OptionData_RAW)
    
    date = OptionData_RAW{i,2};
    date = datenum(int2str(date),'yyyymmdd');
    
    CPflag = char(OptionData_RAW{i,6});
    %strike was 1k-multiplied in original data as provided by OptionMetrics
    strike = OptionData_RAW{i,7}/1000;
    %remove zero-bid quotes
    if(OptionData_RAW{i,8} == 0)
        continue
    end
    spot = IndexData(find(IndexData(:,1)==date),2);
    
    midPrice = (OptionData_RAW{i,8}+OptionData_RAW{i,9})/2; %ask/bid mid
    impVol = OptionData_RAW{i,12};
    
    expiry = datenum(int2str(OptionData_RAW{i,4}),'yyyymmdd');
    
    ticker = cell2mat(OptionData_RAW{i,3});
    
    volume = OptionData_RAW{i,10};
    bidaskSpread = OptionData_RAW{i,9}-OptionData_RAW{i,8};
    OI = OptionData_RAW{i,11};
    last_date = OptionData_RAW{i,5};
    if(isnan(last_date))
        last_date = '19000101';
    end
    last_date = datenum(int2str(last_date),'yyyymmdd');
    
    OptionData_Clean = [OptionData_Clean;[date,CPflag == 'C',isempty(strfind(ticker,'SPXW')),midPrice,expiry,strike,spot,impVol,volume,bidaskSpread,OI,last_date]]; 
    
    i
end

OptionData_Clean = OptionData_Clean(13:end,:);

data = OptionData_Clean;

%read the discounting factor and align to each option quote
treasuryData = dlmread('US 3M Treasury 1996-2018 WRDS Fed.csv',',',1,0);
treasuryDates = treasuryData(:,1);

alignedData = [data -1*ones(size(data,1),1)];
lastRow = 1;

for j= 1:size(alignedData,1)
    for i = lastRow:size(treasuryData,1)
        
        Tdate = datenum(int2str(treasuryData(i,1)),'yyyymmdd');
        
        if(Tdate==alignedData(j,1))
            alignedData(j,13) = treasuryData(i,2) / 100;
            lastRow = i;
            break
        end
    end
    
    j
end

save('SPX Clean 1996-2017 BT Full Maturity','vol30dATM','OptionData_Clean','IndexData_Clean','alignedData')