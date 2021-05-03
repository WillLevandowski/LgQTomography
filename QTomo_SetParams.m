frequencies=[0.75;1.5;3;6;12];
all_at_once=1;
 
minlat = 22.5; 
maxlat = 50; 
maxlong = -79; 
minlong = -125;
minv=4; %%% events must be recorded on minv stations; stations must record minv events
mindist=110;
maxdist=1600;

gamma=0.5;
B = 3.5; %% shear velocity

smooth_coefficient=40; %% a weight. Chosen by Tikhonov regularization
smooth_f_exponent=0.8; %% smoothing matrix :: smooth_coeff * f^exponent

%%% bounds of the Q model
leftlong=-125;%minlong;
rightlong=-79;%maxlong;
toplat=50;
botlat=22.5;
    
%%% node spacing / cell width in the two inversion stages
%%% add an irrational number to one of the spacings to ensure that points
%%% will not coincide, which avoids one of many stability-well pitfall
    space1=100+pi;
    space2=50;
