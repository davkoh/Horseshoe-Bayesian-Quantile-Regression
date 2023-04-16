function [output] = datareshape(data,datasize)
if rem(datasize,100)>0
    disp('Please choose a datasize divisible by 100')
    return
end
if datasize<100
    disp('Please choose a datasize larger than 100');
    return
end
if rem(datasize/100,2)==0 %Since each MC design is 200 observations and we want to use same 100 observations to check forecast error we have to check if the datesize/100 is odd or even.
    output=data(1:100,:,1:100);
    iternum=(datasize/100)/2;
    for i=1:(iternum)
        start=i*100+1;
        ending=(i+1)*100;
        if i==iternum %if datesize/100 is even then we take 100 observations for last iteration
            temp=data(1:100,:,start:ending);
        else
            temp=data(:,:,start:ending);
        end
        output=cat(1,output,temp);
    end
else %if datesize/100 is odd we don't have to handle last iteration as 100 is taken from first 100 designs
    output=data(1:100,:,1:100);
    iternum=floor((datasize/100)/2); %we round down for odd datesize/100
    for i=1:(iternum)
        start=i*100+1;
        ending=(i+1)*100;
        temp=data(:,:,start:ending);
        output=cat(1,output,temp);
    end
end
temp=data(101:200,:,1:100); %taking the data we will run the forecast evaluation on
output=cat(1,output,temp); %appending to the output data
end

