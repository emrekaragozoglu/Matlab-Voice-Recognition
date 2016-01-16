fs=16000;
num=zeros(24000,1);
numarray=zeros(24000,10);

numMFCC=zeros(650,10);
MFCCarray=zeros(650,10);
cc=1;
for i=1:10 %%pronounced numbers (0,1,2,....,9)
for column=1:5 %%number of times each digit is pronounced
    cc
    num = wavrecord(1.5*fs,fs,'double');
    clc
    cc=cc+1;
    
    for j=1:24000
        numarray(j,i)=num(j);
    end
    frame_duration= 0.025;
    frame_len = fs*frame_duration;
    Nstart=0;
    N=10000;
    for n=1:24000
        if(num(n)>0.04)
            Nstart=n;
            break
        end     
    end
    
    numcrop= num(Nstart+1:(Nstart+10000));%25frames
    num_frames = floor (N/frame_len);

    h=floor(hamming(frame_len)*128);
    h(h==128)=127;

    fmel= linspace(0,509,28);
    forg=700*((exp(fmel/1125))-1);
    forgfloor=floor(forg);

    dctlogwindowedpower = zeros(1,26);
    count=1;
    for a=1:num_frames
        frame = numcrop((a-1)*frame_len+1 :frame_len*a);

        framei = frame.*h;
        fftframei= floor((fft(framei,frame_len))*128);
        powerframei=((abs(fftframei)).^2)/(frame_len);

        for k= 1:1:26

            if(k<27)
                tri = floor(triang(forgfloor(k+2)-forgfloor(k))*128);
            end
            if(k==1)
              win=[tri;zeros(400-length(tri),1)];
            elseif(k==26)
              win=[zeros(forgfloor(k),1);tri];  
            else
              if(k<27)
                win=[zeros(forgfloor(k),1);tri;zeros(400-forgfloor(k+2),1)]; 
              end
            end    

            windowedpower=powerframei.*win;
            sumwindowedpower=0;    
            if(k==1)
              for s=1:1:forgfloor(k+2)
                  sumwindowedpower=sumwindowedpower+windowedpower(s);
              end   
                logwindowedpower=log(sumwindowedpower);
            elseif(k<27)
              for s=forgfloor(k)+1:1:(forgfloor(k+2))
                    sumwindowedpower=sumwindowedpower+windowedpower(s);
              end    
                   logwindowedpower=log(sumwindowedpower);
            end
              dctlogwindowedpower=dct(logwindowedpower);
              numMFCC(count,column)=dctlogwindowedpower;
              count = count +1;

        end
    end
    
              sum2=0;
              average=0;
              %%%averaging

              for row=1:650
                  for col=1:10
                    sum2=sum2+numMFCC(row,col);  
                  end
                  average=sum2/5;
                  MFCCarray(row,i)=average;
                  sum2=0;
              end
end
    cc=1;
end

