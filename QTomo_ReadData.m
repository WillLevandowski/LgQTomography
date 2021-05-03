%%% reads in data from text files  
%%% culls data to chosen lat/lon bounds and distance range
%%% culls events and stations with too few (<"minv") data

%%% Next, accounts for geometric spreading 
%%% Here, also computes log(amplitude) residuals realtive to the best-fit
%%% constant-Q, constant-S, constant-R model


%%% This input section is contrived and will probably need to be re-written
%%% unless data is formatted exactly as in "Q_data_0.75.txt" and others.

   if f==0.75
       [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp] = textread('Qdata_0.75.txt','%s %f %f %s %f %f %f %f %f %f %f');
   end
   if f==1.5
       [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp] = textread('Qdata_1.5.txt','%s %f %f %s %f %f %f %f %f %f %f');
   end
   if f==3
       [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp] = textread('Qdata_3.txt','%s %f %f %s %f %f %f %f %f %f %f');
   end
   if f==6
       [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp] = textread('Qdata_6.txt','%s %f %f %s %f %f %f %f %f %f %f');
   end
   if f==12
       [event,evlat,evlong,sta,stlat,stlong,dist,time,lfq,hfq,amp] = textread('Qdata_12.txt','%s %f %f %s %f %f %f %f %f %f %f');
   end

   
   %%% Now is a good time to enumerate any troublesome stations and events
   %%% Setting dist>distmax will ensure they are discarded.
  %%%% this station had some problems, so it's excluded.
   foo=find(stlong==-104.7060);dist(foo)=1e10;
%    foo=find(stlong<-111 & dist>1300);dist(foo)=1e10; 
    
   %%% chop arrivals from stations and events outside of chosen area or
   %%% chosen distance range. (see set_params.m)
   loc = find(evlat>=minlat & evlat<=maxlat & stlat >= minlat & stlat <= maxlat ...
            & evlong >= minlong & evlong <= maxlong & stlong >= minlong & stlong <= maxlong ...
         & dist>mindist & dist<maxdist );
    event = event(loc);
    sta = sta(loc);
    dist = dist(loc);
    amp = amp(loc);
    stlong=stlong(loc);
    evlong=evlong(loc);
    stlat=stlat(loc);
    evlat=evlat(loc);
    
%%% plot the amplitudes. Optional.
%      figure(1000+ceil(f));semilogy(dist,amp,'k.');hold on;title([ num2str(f) ' Hz amplitudes']);set(gca,'YTick',[1e-6 1]); yticks([1e-6 1])
   
    stations = unique(sta); 
    sources = unique(event);
    nso = length(sources);
    nsta = length(stations);  
    namp = length(amp);                  
    row = nso*nsta;
    col = nso + nsta;      % Determine number of elements in G matix
    
    sonum=zeros(namp,nso);
    stanum=zeros(namp,nsta);
    for i = 1:nso;  sonum(:,i) = strcmp(sources(i),event);  end  % Create matrix for Sources   
    for i = 1:nsta; stanum(:,i) = strcmp(stations(i),sta); end     % Create matrix for Receivers
    
    Gstasum = sum(stanum); 
    Gsta = find(Gstasum>=minv);
    goodsta = stations(Gsta);
    tsta = ismember(sta,goodsta);
    
    %%% check that stations have enough arrivals and events have enough
    %%% recordings, at least "minv" (see set_params.m)
    Gsosum = sum(sonum); 
    Gsor = find(Gsosum>=minv);
    goodsor = sources(Gsor);
    tsor = ismember(event,goodsor);
    newloc = find(tsta == 1 & tsor == 1);
    
    %%% cull data 
    event2 = event(newloc);
    sta2 = sta(newloc);
    dist2 = dist(newloc);
    amp2 = amp(newloc);
    stlong2=stlong(newloc);
    evlong2=evlong(newloc);
    stlat2=stlat(newloc);
    evlat2=evlat(newloc);
    
    stations2 = unique(sta2);
    sources2 = unique(event2);
    nso2 = length(sources2);
    nsta2 = length(stations2);
    namp2 = length(amp2);                       
    
    row2 = nso2*nsta2;col2 = nso2 + nsta2;    % Determine number of elements in G matix
    
    sonum2=zeros(namp2,nso2);
    stanum2=zeros(namp2,nsta2);
    for i = 1:nso2;sonum2(:,i) = strcmp(sources2(i),event2);end
    for i = 1:nsta2;stanum2(:,i) = strcmp(stations2(i),sta2);end
    
    
    %%% "amp_" is log{amplitude corrected for geometric spreading}
    amp_=log(amp2) + gamma.*log(dist2*1000);    
    p=polyfit(dist2,amp_,1); %%% p(1) is slope, proportional to Q. 
    %%% p(2) is intercept. Physically, p(2) is the average 
    %%% (source term + receiver term). Since average R~0, p(2)~average S.
    %%% for the data at this frequency. 
    
    
   %%% Optional. Check spreading-corrected amplitude residuals to identify problematic data.
%     pred=polyval(p,dist2);
%     res=amp_-pred;
%     for example, foo=find(abs(res)>quantile(abs(res),0.99)); 
%     
    
%%% generate list of lat/lon of stations with data at this frequency
    [sta_locations,~]=unique([stlong2 stlat2],'rows');

    %%% Compute the "signal", which is the vector of amplitude variations
    %%% with respect to a constant Q, constant S, constant R==0 reference
    %%% model. Results don't depend on this step as long as the smoothing
    %%% parameter, node spacing, etc., are tuned appropriately for a stable
    %%% inversion. Gutting the background effects can set up a more stable
    %%% inversion without needing extra damping, a logarithmic barrier to
    %%% prevent negsative 1/Q (or similar). Those additional G matrix
    %%% components can slow down the inversion too much or push it over
    %%% default memory allocation.
    
%%% "signal" is the remaining log(amplitude) after accouting for geometric
%%% spreading, average dataset-wide attenuation, and the average intercept.
    background_attenuation= p(1)*B/(-pi*f); %%% best-fitting 1/Q
    background_source=p(2); %%% reference source term at this frequency
     signal = amp_ - background_source + pi*f*dist2*background_attenuation/B;
  