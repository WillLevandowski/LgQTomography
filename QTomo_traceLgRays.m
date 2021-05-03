%%% create the model geometry.
%%% make a cartesian model from lat/lon.
%%% it will be centered (east-west) on the center of the study

%%% distances in km

%%% this script is a mess and is poorly commented. It was written five
%%% years ago and basically not altered since, so some of the logic is
%%% hazy. 

centerlong=leftlong/2+rightlong/2;
centerlat=toplat/2+botlat/2;
xmax=(rightlong-centerlong)*111.12*cosd(botlat);
xmin=(leftlong-centerlong)*111.12*cosd(botlat);
ymax=(toplat-centerlat)*111.12;
ymin=(botlat-centerlat)*111.12;

All_Events=[evlong2 evlat2 stlong2 stlat2 ];

ev_Y=(All_Events(:,2)-centerlat)*111.12;
ev_X=(All_Events(:,1)-centerlong).*(111.12*cosd(All_Events(:,2)));
st_Y=(All_Events(:,4)-centerlat)*111.12;
st_X=(All_Events(:,3)-centerlong).*(111.12*cosd(All_Events(:,4)));

% total_pathlengths=sqrt( (ev_X-st_X).^2 + (ev_Y-st_Y).^2);

xs=(xmin:spacing:xmax+spacing/2);
ys=(ymin:spacing:ymax+spacing/2);
[x,y]=meshgrid(xs,ys);
X=reshape(x,length(xs)*length(ys),1);
Y=reshape(y,length(X),1);
Lat=centerlat+Y/111.12;
Lon=centerlong+X./(111.12*cosd(Lat));

%%% trim up the edges, since inversion time scales ~with length(X)^2.
%%% Also, this step is a chance to get rid of unresolvable nodes
okay=find(Lon>=leftlong-0.5*spacing/111.12 ...
    & Lon<=rightlong+0.5*spacing/111.12 ...
    & Lat<=toplat+0.5*spacing/111.12 ...
    & Lat>=botlat-0.5*spacing/111.12...
   ); 
Lon=Lon(okay);
Lat=Lat(okay);
X=X(okay);
Y=Y(okay);


%%%% with the model geometry set, trace the raypaths

rad=spacing*sqrt(2); %% limit the area searched for nodes
hits=0*X; %%% hit count in each cell (will be weighted)
hit=zeros(length(X),length(All_Events)); %% does arrival i sense node j?
pathlength=zeros(length(X),length(All_Events)); %% similar to hit

%%% McNamara's approach (after Benz et al., 1997?) uses the RMS spectral
%%% amplitude of 'displacement'. 'Displacement' is the time series over the
%%% window attendant to Lg. Window begins with 3.6 km/s arrival and ends
%%% with 3.0 km/s. For example, if hypocentral distance = 1080 km, window
%%% runs from 300 s after origin to 360 s after origin.
% VLg_min=3.0; %%% velocity for end time window
% VLg_max=3.6; %%% velocity for start of time window


for j=1:length(ev_X) %%% foreach ray
    
    
    %%% this step is specific to the McNamara time-window approach
    %%% time window = origin + dist/3.6 until origin + dist/3.0
%     dtmax=dist2(j)/VLg_min-dist2(j)/VLg_max; 
%     dlmax=dtmax/2*B;
%     lmax=dist2(j)/2+dlmax;
%     dxmax=sqrt(lmax^2 - (dist2(j)/2)^2);%%% max distance away from ray that gets to receiver within the time duration
    
    %%% I haven't checked, but it might be possible to hack out 
    %%% finite-width kernels by setting dxmax to something much smaller 
    %%% than node spacing.
    
    %%% This section is ugly and almost certainly not optimized. But it
    %%% works, so I'll keep it for now
    
    %%% ray's start and end coordinates
    x1=ev_X(j);
    x2=st_X(j);
    y1=ev_Y(j);
    y2=st_Y(j);
    
    for i=1:length(X)
        
        %%%% check each node/cell
        x0=X(i);
        y0=Y(i);
        d=abs( (y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2);
        
        
        %%% foreach node/cell near to or traversed by the ray, 
        %%% calculate the radial sensitivity kernel about the cell. At any
        %%% point along the raypath, the kernel is annu
        %%% the sensit
        if d<rad && ...
           ((x0<x1+spacing/2 && x0>x2-spacing/2) ||  (x0>x1-spacing/2 && x0<x2+spacing/2) ) && ...
           ((y0<y1+spacing/2 && y0>y2-spacing/2) ||  (y0>y1-spacing/2 && y0<y2+spacing/2) )

             % What portion of the path has been covered?
             dx=x2-x1; 
             dxleft=(x0-spacing/2)-x1; 
             dxright=(x0+spacing/2)-x1;
             x_left=x1+dxleft;
             x_right=x1+dxright;
             portion_right=dxright/dx;
             portion_left=dxleft/dx;

             dy=y2-y1; 
             dybot=(y0-spacing/2)-y1;
             dytop=(y0+spacing/2)-y1;
             y_bot=y1+dybot;
             y_top=y1+dytop;
             portion_bot=dybot/dy;
             portion_top=dytop/dy;

             y_right=dy*portion_right+y1;
             y_left=dy*portion_left+y1;
             x_top=dx*portion_top+x1;
             x_bot=dx*portion_bot+x1;

             right=[x_right y_right];
             left=[x_left y_left];
             top=[x_top y_top];
             bot=[x_bot y_bot];
             corners=[top;left;bot;right;];
             portions=[portion_top;portion_left;portion_bot;portion_right;];


             %%% where does the ray enter and leave ("pierce") the cell?

             piercing=find( abs(corners(:,1)-x0)<=spacing/2+0.001 & abs(corners(:,2)-y0)<=spacing/2+0.001);
             if abs(x0-x2)<spacing/2 && abs(y0-y2)<spacing/2
                 piercing= find( abs(corners(:,1)-x0)<=spacing/2+0.001 & abs(corners(:,2)-y0)<=spacing/2+0.001 & portions>0); 
                 corners=[corners(piercing,:); x2 y2];
                 piercing=[1;2];
             end
             if abs(x0-x1)<spacing/2 && abs(y0-y1)<spacing/2 
                 corners=[top;left;bot;right;];
                 piercing=find( abs(corners(:,1)-x0)<=spacing/2+0.001 & abs(corners(:,2)-y0)<=spacing/2+0.001 & portions>0); 
                 corners=[corners(piercing,:); x1 y1];
                 piercing=[1;2];
             end
             while length(piercing)>2
                 rand_corner=ceil(rand*length(corners));
                 repeat=find(abs(corners(:,1)-corners(rand_corner,1))<1 & abs(corners(:,2)-corners(rand_corner,2))<1);
                 if length(repeat)>1
                     piercing(repeat)=rand_corner;
                     piercing=unique(piercing);
                 end
             end
             %%% the entry/exit points have been identified
             if length(piercing)==2
                    in_cell=sqrt( (corners(piercing(1),1)-corners(piercing(2),1))^2 + (corners(piercing(1),2)-corners(piercing(2),2))^2);

                    pierce_points1=[corners(piercing(1),1) corners(piercing(1),2)];
                    pierce_points2=[corners(piercing(2),1) corners(piercing(2),2)];
                    midpoint=pierce_points1/2+pierce_points2/2;

%                     d1=sqrt( (midpoint(1)-ev_X(j)).^2 + (midpoint(2)-ev_Y(j)).^2 );
%                     d2=sqrt( (midpoint(1)-st_X(j)).^2 + (midpoint(2)-st_Y(j)).^2 );
% 
%                     d_=2*min([d1 d2])/dist2(j); %%% proportion of raypath
% 
%                     dx=dxmax*d_; %%% max off-raypath distance for a forward-scatterer to contribute energy to the window
% 
%                     p_total=  sqrt( (X-st_X(j)).^2 +  (Y-st_Y(j)).^2 ) + ...
%                                 sqrt( (X-ev_X(j)).^2 +  (Y-ev_Y(j)).^2 );
%                     d=sqrt( (X-midpoint(1)).^2+ (Y-midpoint(2)).^2);
% 
%                     %%% dp is extra pathlength required for
%                     %%% source-scatter-receiver ray relative to direct
%                     %%% source-receiver path
%                     dp=p_total-dist2(j);
%                     dt=dp/B; % time delay relative to direct ray
%                         
%                     %%% all the nodes that would contribute
%                     %%% forward-scattered energy to the time window
%                     near=find(dt<dtmax & d<dxmax);
%                     w=1-(dt(near)/dtmax);
%                     w=w/sum(w);
%                     
%                     if length(near)<4
                          near=find(abs(X-midpoint(1))<spacing & abs(Y-midpoint(2))<spacing );
                          d=sqrt( (X(near)-midpoint(1)).^2+ (Y(near)-midpoint(2)).^2);
                          w=1./d;w=w/sum(w);
%                     end

                pathlength(near,j)= pathlength(near,j)+in_cell*w;
                hits(near)=hits(near)+w;
                hit(i,j)=1; 
             end
        end
    end
end

%%% normalized pathlengths, 'path_proportion'-->"path_prop"
path_prop=0*pathlength;
for i=1:length(ev_X)
    path_prop(:,i)=pathlength(:,i)/sum(pathlength(:,i));
end
path_prop=path_prop';

%   last = form*dist2; %%% This quantity may not be used.
path_props=0*path_prop;
for i=1:length(amp2)
    path_props(i,:)=path_prop(i,:)*form*dist2(i);
end

%%%%% some metrics
% mean_path=0*X;
% median_path=0*X;
% smr_path=0*X;
% for i=1:length(X)
%     paths_through=find(hit(i,:));
%     dist_through=dist2(paths_through);
%     mean_path(i)=mean(dist_through);
%     median_path(i)=median(dist_through);
%     smr_path(i)=mean(sqrt(dist_through))^2;
% end  
% 
% dlmwrite('mean_path', [Lon Lat mean_path], '\t')
% dlmwrite('median_path', [Lon Lat median_path], '\t')
% dlmwrite('smr_path', [Lon Lat smr_path], '\t')
