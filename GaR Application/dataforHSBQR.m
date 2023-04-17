% This code is an edited version of McCraken's "fredfactors.m". In
% particular it is the first part of the fredfactors code edited for the
% HSBQR

% PART 1: LOAD AND LABEL DATA
csv_in='current.csv';
% Load data from CSV file
dum=importdata(csv_in,',');

% Variable names
series=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(2,:); %edited because tcode is 3rd row of current.csv
factorinc=dum.data(1,:); %will we use smaller database?

% Raw data
rawdata=dum.data(3:end,:); %edited because data is from 4th row of current.csv

% Month/year of final observation
final_datevec=datevec(dum.textdata(end,1));
final_month=final_datevec(2);
final_year=final_datevec(1);

% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
dates = (1959+3/12:3/12:final_year+final_month/12)'; %edited for quarterly data

% T = number of months in sample
T=size(dates,1);
rawdata=rawdata(1:T,:);

% =========================================================================
% PART 2: PROCESS DATA

% Transform raw data to be stationary using auxiliary function
% prepare_missing()
ytraw=prepare_missing(rawdata,tcode);
ytraw=ytraw(2:T,:); %first difference creates NA, so we remove those entries
dates=dates(2:T,:);

widedata=[]; %We chose a data structure that starts in 1970Q1 as that gives a large selection of variables
for i=2:size(ytraw,2)
    temp=isnan(ytraw(44,i));
    if temp==0
        datatemp=ytraw(:,i);
        widedata=[widedata,datatemp];
    end
end

lastsum=isnan(sum(widedata(end,:))); %Check if last row has NaN's. If they do we omit them from analysis
if lastsum==1
    widedata=widedata(44:(end-1),:);
    ywide=ytraw(44:(end-1),1);
else
    widedata=widedata(44:end,:);
    ywide=ytraw(44:end,1);
end

% Standardising Data
wideX=[];
for i=1:size(widedata,2)
    m=mean(widedata(:,i));
    st=std(widedata(:,i));
    x_a=(widedata(:,i)-m)./st;
    wideX=[wideX,x_a];
end

wideX=[ones(size(wideX,1),1), wideX];
wideDB=[ywide,wideX];

clearvars -except wideDB
save('GaR application\processed\GaRdata.mat','wideDB'); %saving as a .mat file
%writematrix(wideDB,'Empirical application\processed\wideDB.xlsx')
xlswrite('GaR application\processed\wideDB.xlsx',wideDB)