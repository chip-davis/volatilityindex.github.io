clear

%====Load Data====
load('SPX Clean 1996-2017 BT Full Maturity');
load('SPX ATM Vol30d')

%====Initialize Variables====
dates = unique(alignedData(:,1));
tmp = zeros(length(dates),6);
tmp2 = zeros(length(dates),40);
tmp3=[];
phi_storage = {};
result = zeros(length(dates),11);
weekday_ts = weekday(dates);
refresh = true;
lastDay = 1;
putCount = 0;
callCount = 0;



%Before this day alpha is estimated weekly; afterwards daily
modeSwitch = 4278;

%Estimate alpha^pm
for i = 1:length(dates) 
    
    
    
    %Find the data of day i and all available expiries
    date = dates(i);
    todayData_all = alignedData(alignedData(:,1)==date,:);
    expiry = todayData_all(:,5);
    expiryList = unique(expiry);
    
    %Re-Initialize temporary variables in the beginning of a new week
    if(refresh == true)
        alpha_t_minus_est = zeros(1,1);
        putDataSummary = zeros(1,size(alignedData,2));
        alpha_t_plus_est = zeros(1,1);
        callDataSummary = zeros(1,size(alignedData,2));
        putCount = 0;
        callCount = 0;
        refresh = false;
    end
    
    phi_t_minus_est = zeros(1,3);
    phi_t_plus_est = zeros(1,3);
        
    %Loop through expiries
    %Store alpha estimation points in alpha_t_plus_est and alpha_t_minus_est
    for j = 1:length(expiryList)
        
        todayData = todayData_all(todayData_all(:,5)==expiryList(j),:);
        
        callData = todayData(todayData(:,2)==1,:);
        callData = sortrows(callData,6);
        putData = todayData(todayData(:,2)==0,:);
        putData = sortrows(putData,6);
        
        %Proceed only if there are at list 2 valid data points
        if(size(callData,1)<=1 || size(putData,1)<=1)
            continue;
        end
        
        %Proceed only if maturity is within [6,31] trading days
        maturityTradingDays = length(busdays(callData(1,1),callData(1,5)));
        
        if(maturityTradingDays<6 || maturityTradingDays>31) 
            continue;
        end
        
        
        %Interpolate ATM implied volatility
        %Select only the deep OTM quotes
        tau_C = (length(busdays(callData(:,1),callData(:,5)))/252);% * ones(size(callData,1),1); 
        tau_P = (length(busdays(putData(:,1),putData(:,5)))/252);% * ones(size(putData,1),1);
        
        fwdPrice = exp(-callData(1,end)* tau_C(1))*callData(1,7);
        callData = callData(callData(:,6)>=fwdPrice,:);
        minCall = callData(callData(:,6)==min(callData(:,6)),:);
        putData = putData(putData(:,6)<=fwdPrice,:);
        maxPut = putData(putData(:,6)==max(putData(:,6)),:);
        
        %Proceed only if there are at list 2 valid data points
        if(size(callData,1)<=1 || size(putData,1)<=1)
            continue;
        end
        
        impVol = ((fwdPrice-maxPut(1,6))*minCall(1,8)+(minCall(1,6)-fwdPrice)*maxPut(1,8))/(minCall(1,6)-maxPut(1,6));
                
        tmp2(i,j) = impVol;
        tmp3= [tmp3; [impVol ((maxPut(1,6))*minCall(1,8)+(minCall(1,6))*maxPut(1,8))/(minCall(1,6)+maxPut(1,6)) (fwdPrice-maxPut(1,6))/(minCall(1,6)-maxPut(1,6)) maxPut(1,6)/(minCall(1,6)+maxPut(1,6)) ]];
        
        
        %====Begin: remove quotes which violate no-arbitrage====
        % See documentation for no-arbitrage rule details
        tempPx = putData(end,4);
        omitPut = zeros(size(putData,1),1);
        for l = (size(putData,1)-1):(-1):1
            px = putData(l,4);
            if(px - tempPx >= -1e-10)
                omitPut(l,1) = 1;
                continue;
            end
            tempPx = px;
        end
        putData = putData(omitPut==0,:);

        tempPx = callData(1,4);
        omitCall = zeros(size(callData,1),1);
        for l = 2:size(callData,1)
            px = callData(l,4);
            if(px - tempPx >= -1e-10)
                omitCall(l,1) = 1;
                continue;
            end
            tempPx = px;
        end
        callData = callData(omitCall==0,:);        
        %====End: remove quotes which violate no-arbitrage====
        
        mC = log(callData(:,6)./fwdPrice);
        logMoneynessC = mC ./ (impVol * sqrt(tau_C));
        mP = log(putData(:,6)./fwdPrice);
        logMoneynessP = mP ./ (impVol * sqrt(tau_P));
        
        callData = callData(logMoneynessC>=1,:);
        mC = mC(logMoneynessC>=1);
        logMoneynessC = logMoneynessC(logMoneynessC>=1);
        putData = putData(logMoneynessP<=-2.5,:);
        mP = mP(logMoneynessP<=-2.5);
        logMoneynessP = logMoneynessP(logMoneynessP<=-2.5);
        
        Nt_C=size(callData,1);
        Nt_P=size(putData,1);
        
        %Store price-moneyne slope of each option pair in alpha_t_plus_est
        %and alpha_t_minus_est
        if(size(callData,1)>1)
            alpha_t_plus_est = [alpha_t_plus_est ; log(callData(1:(Nt_C-1),4)./callData(2:Nt_C,4))./(mC(1:(Nt_C-1),1)-mC(2:Nt_C,1))];
            callCount = callCount + size(callData,1);
            callDataSummary = [callDataSummary; callData];
        end
        if(size(putData,1)>1)
            alpha_t_minus_est =[alpha_t_minus_est ; log(putData(1:(Nt_P-1),4)./putData(2:Nt_P,4))./(mP(1:(Nt_P-1),1)-mP(2:Nt_P,1))];
            putCount = putCount + size(putData,1);
            putDataSummary = [putDataSummary; putData];
        end
        
        if(size(callData,1)>=3)
            tau_C = (length(busdays(callData(:,1),callData(:,5)))/252) * ones(size(callData,1),1); 
            phi_t_plus_est = [phi_t_plus_est;[log(exp(callData(:,end).*tau_C).*callData(:,4)./(tau_C .* exp((callData(:,end)).* tau_C).*callData(:,7))) mC tau_C]];
        end
        if(size(putData,1)>=3)
            tau_P = (length(busdays(putData(:,1),putData(:,5)))/252) * ones(size(putData,1),1);
            phi_t_minus_est = [phi_t_minus_est; [log(exp(putData(:,end).*tau_P).*putData(:,4)./(tau_P .* exp((putData(:,end)).* tau_P).*putData(:,7))) mP tau_P]];
        end
    end
    
    %Store phi^pm factors for later calculation
    phi_storage{i,1} = phi_t_minus_est(2:end,:);
    phi_storage{i,2} = phi_t_plus_est(2:end,:);
    
    %Estimate alpha at the end of week (on or before mode switch) or each
    %day (after mode switch); Estimate only if there are at least three
    %valid option pairs; Use previous alpha estimation if otherwise
    if(i > 1 && i < modeSwitch && weekday_ts(i) > weekday_ts(i+1))
        alpha_t_plus_est = alpha_t_plus_est(2:end);
        alpha_t_minus_est = alpha_t_minus_est(2:end);
        alpha_t_plus = 1-median(alpha_t_plus_est);
        alpha_t_minus = -1+median(alpha_t_minus_est);
                
        
        if(lastDay > 1 && (isnan(alpha_t_plus)||size(alpha_t_plus_est,1)<3))
            alpha_t_plus = result(lastDay-1,5);
        end
        
        if(lastDay > 1 && (isnan(alpha_t_minus)||size(alpha_t_minus_est,1)<3))
            alpha_t_minus = result(lastDay-1,4);
        end
        
        result(lastDay:(i-1),:)= [dates(lastDay:(i-1)) putCount*ones(i-lastDay,1) callCount*ones(i-lastDay,1) alpha_t_minus*ones(i-lastDay,1) alpha_t_plus*ones(i-lastDay,1) zeros(i-lastDay,6)];
        
        lastDay = i;
        refresh = true;
    elseif(i == modeSwitch)
        alpha_t_plus_est = alpha_t_plus_est(2:end);
        alpha_t_plus = 1-median(alpha_t_plus_est);
        alpha_t_minus_est = alpha_t_minus_est(2:end);
        alpha_t_minus = -1+median(alpha_t_minus_est);
        
        if(lastDay > 1 && (isnan(alpha_t_plus)||size(alpha_t_plus_est,1)<3))
            alpha_t_plus = result(lastDay-1,5);
        end
        
        if(lastDay > 1 && (isnan(alpha_t_minus)||size(alpha_t_minus_est,1)<3))
            alpha_t_minus = result(lastDay-1,4);
        end
        
        result(lastDay:i,:)= [dates(lastDay:i) putCount*ones(i-lastDay+1,1) callCount*ones(i-lastDay+1,1) alpha_t_minus*ones(i-lastDay+1,1) alpha_t_plus*ones(i-lastDay+1,1) zeros(i-lastDay+1,6)];
        refresh = true;
    elseif(i >= modeSwitch)
        alpha_t_plus_est = alpha_t_plus_est(2:end);
        alpha_t_plus = 1-median(alpha_t_plus_est);
        alpha_t_minus_est = alpha_t_minus_est(2:end);
        alpha_t_minus = -1+median(alpha_t_minus_est);

        if(isnan(alpha_t_plus)||size(alpha_t_plus_est,1)<3)
            alpha_t_plus = result(i-1,5);
            tmp(i,6)=1;
        end
        
        if(isnan(alpha_t_minus)||size(alpha_t_minus_est,1)<3)
            alpha_t_minus = result(i-1,4);
            tmp(i,5)=1;
        end
        
        result(i,1:5)= [dates(i) putCount callCount alpha_t_minus alpha_t_plus];
        
        refresh = true;
    end
end

%Estimate phi^pm, LJV and other tail risk measures

for i = 1:length(dates)
       
    date = dates(i);
    todayData_all = alignedData(alignedData(:,1)==date,:);
    expiry = todayData_all(:,5);
    expiryList = unique(expiry);
    
    phi_t_minus_est = phi_storage{i,1};
    phi_t_plus_est = phi_storage{i,2};
       
    alpha_t_minus = result(i,4);
    alpha_t_plus = result(i,5);
    
    %Estimate phi^pm with estimated alpha^pm series 
    %and previously stored phi factors
    phi_t_minus = exp((log(alpha_t_minus+1)+log(alpha_t_minus)) +  median(phi_t_minus_est(:,1)-(1+alpha_t_minus)*phi_t_minus_est(:,2)));
    phi_t_plus = exp((log(alpha_t_plus-1)+log(alpha_t_plus)) + median(phi_t_plus_est(:,1)-(1-alpha_t_plus)*phi_t_plus_est(:,2)));
    
    if(size(phi_t_plus_est,1)<=1)
        phi_t_plus = NaN;
    end
    
    if(size(phi_t_minus_est,1)<=1)
        phi_t_minus = NaN;
    end
        
    LJI = NaN;
    LJV = NaN;
    RJI = NaN;
    RJV = NaN;
    
    tP = NaN;
    tC = NaN;
    
    impVol30d = vol30dATM{i,5};
    
    %threshold of large jump: 10 times 30d ATM implied vol
    k_t_P = 10 * impVol30d * sqrt(5/252);
    k_t_C = 10 * impVol30d * sqrt(5/252);
        
    
    
    if(~isnan(phi_t_minus))
        LJI = phi_t_minus * exp(-k_t_P* alpha_t_minus)/alpha_t_minus;
        LJV = phi_t_minus * exp(-alpha_t_minus*k_t_P)*(alpha_t_minus*k_t_P*(alpha_t_minus*k_t_P+2)+2)/(alpha_t_minus^3);    
    end
    if(~isnan(phi_t_plus))
        RJI = phi_t_plus * exp(-k_t_C * alpha_t_plus)/alpha_t_plus;
        RJV = phi_t_plus * exp(-alpha_t_plus*k_t_C)*(alpha_t_plus*k_t_C*(alpha_t_plus*k_t_C+2)+2)/(alpha_t_plus^3);
    end
               
        
    result(i,6:11)= [phi_t_minus phi_t_plus LJI RJI LJV RJV]';
end

calendar = datetime(result(:,1),'convertfrom','datenum');

leftDensity = (5/252).*result(:,6).*exp(-5.*result(:,4).*vol30dATM.impl_volatility.*sqrt(5/252))./result(:,4);
leftDensityFixed = (5/252).*result(:,6).*exp(-0.1.*result(:,4))./result(:,4);

maLength = 21;
leftDensityMA = zeros(length(calendar),1);
leftDensityFixedMA = zeros(length(calendar),1); 
LJVMA = zeros(length(calendar),1); 

for j = maLength:length(calendar)
   leftDensityMA(j) = nanmean(leftDensity((j-maLength+1):j));

   LJVMA(j) = nanmean(result((j-maLength+1):j,10));


   if(sum(isnan(leftDensityFixed((j-maLength+1):j)))<=maLength/2)
       leftDensityFixedMA(j) = nanmean(leftDensityFixed((j-maLength+1):j));
   else
       leftDensityFixedMA(j) = nan;
   end

end


save('SPX Tail Risk Measures','calendar','result','LJVMA','leftDensityFixedMA');