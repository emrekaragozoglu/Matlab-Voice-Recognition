    fs=16000;
    numtest=zeros(24000,1);
    numtestMFCC=zeros(650,1);
    numtest = wavrecord(1.5*fs,fs,'double');

    frame_duration= 0.025;
    frame_len = fs*frame_duration;
    Nstart=0;
    N=10000;
    for n=1:24000
        if(numtest(n)>0.04)
            Nstart=n;
            break
        end     
    end
    
    numtestcrop= numtest(Nstart+1:(Nstart+10000));%25frames
    wavplay(numtestcrop,fs);
    
    num_frames = floor (N/frame_len);

    h=floor(hamming(frame_len)*128);
    h(h==128)=127;

    fmel= linspace(0,509,28);
    forg=700*((exp(fmel/1125))-1);
    forgfloor=floor(forg);

    dctlogwindowedpower = zeros(1,26);
    count=1;
    for a=1:num_frames
        frame = numtestcrop((a-1)*frame_len+1 :frame_len*a);

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
%             plot(win);
%             hold on;

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
              numtestMFCC(count)=dctlogwindowedpower;
              count = count +1;

        end
    end
    
    difference=0;
    square=0;
    sum=0;
    sumnum=zeros(10,1);
    %%comparison
    for c=1:10
        for r=1:650
            difference=MFCCarray(r,c)-numtestMFCC(r,1);
            sum=sum+(difference*difference);
        end
        sumnum(c,1)=sum;
        sum=0;
    end
    
[value,index]=min(sumnum);
if(index==10)
    number=0;
else
    number=index;
end
    number
