
%%%% write results

%%%% it is a good idea to write an objective script for color palettes.
%%%% Ex/    quantile(Q,0.01)-->darkred 
%%%%        quantile(Q,0.16)-->red 
%%%%        median(Q)-->white
%%%%        quantile(Q,0.84)-->blue 
%%%%        quantile(Q,0.99)-->darkblue
%%%% Fixing the contours to the quantiles of Q structure facilitates
%%%% comparison of Q across frequencies


%%%% match R and S to station and source ID and to their locations
source_longitude=zeros(nso2,1)+NaN;
source_latitude=zeros(nso2,1)+NaN;
for i=1:nso2
        match=find(strcmp(event2,sources2(i,:)),1,'first');
        source_longitude(i)=evlong2(match);
        source_latitude(i)=evlat2(match);
end

receiver_longitude=zeros(nsta2,1)+NaN;
receiver_latitude=zeros(nsta2,1)+NaN;
for i=1:nsta2
        match=find(strcmp(sta2,stations2(i,:)),1,'first');
        receiver_longitude(i)=stlong2(match);
        receiver_latitude(i)=stlat2(match);
end


%%%%% Remove unconstrained nodes or set them to NaN (or retain),
%%%%%  depending on how the models will be used.
%%%%%  To calculate Q0 and n, remove unconstrained nodes.
Qtemp= Q;
Qtemp(hits==0)=NaN;
  if f==0.75
        dlmwrite(['Inv1_Q_0_75'], [Lon Lat Qtemp], '\t')
        dlmwrite(['Inv1_R_0_75'], [receiver_longitude receiver_latitude dR], '\t')
        dlmwrite(['Inv1_S_0_75'], [source_longitude source_latitude dS+background_source], '\t')
   end
   if f==1.5
        dlmwrite(['Inv1_Q_1_5'], [Lon Lat Qtemp], '\t')
        dlmwrite(['Inv1_R_1_5'], [receiver_longitude receiver_latitude dR], '\t')
        dlmwrite(['Inv1_S_1_5'], [source_longitude source_latitude dS+background_source], '\t')

   end
   if f==3
        dlmwrite(['Inv1_Q_3'], [Lon Lat Qtemp], '\t')
        dlmwrite(['Inv1_R_3'], [receiver_longitude receiver_latitude dR], '\t')
        dlmwrite(['Inv1_S_3'], [source_longitude source_latitude dS+background_source], '\t')
   end
   if f==6
        dlmwrite(['Inv1_Q_6'], [Lon Lat Qtemp], '\t')
        dlmwrite(['Inv1_R_6'], [receiver_longitude receiver_latitude dR], '\t')
        dlmwrite(['Inv1_S_6'], [source_longitude source_latitude dS+background_source], '\t')
   end
   if f==12
        dlmwrite(['Inv1_Q_12'], [Lon Lat Qtemp], '\t')
        dlmwrite(['Inv1_R_12'], [receiver_longitude receiver_latitude dR], '\t')
        dlmwrite(['Inv1_S_12'], [source_longitude source_latitude dS+background_source], '\t')
   end

    
    
    