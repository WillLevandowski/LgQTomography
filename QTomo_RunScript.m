clear all
close all

QTomo_SetParams

for freq_invert=1:length(frequencies)
        %%% set final constants
        f=frequencies(freq_invert);
        form = -pi*f/B; %%% coefficient with 1/Q in exponent 
        spacing=space1;
    
    QTomo_ReadData %%% import data and format
    fprintf(['Inversion at ' num2str(f) ' Hz: \n' ])%...
    QTomo_traceLg %%% build model geometry and trace rays/compute kernels
    QTomo_BuildSmoothingMatrix  %%% build "smooth_A", the smoothing matrix
   tic
    QTomo_Inversion %%% run the inversion
    toc
    QTomo_WriteResults    

end

%%% An example routine to calculate Q0 and eta from 
%%%     Q(f)_node = Q0_node * f ^ (eta_node)

%%% assuming that models at each frequency are at the same locations
Qf=zeros(length(Lon),length(frequencies));
dat=load('Inv1_Q_0_75'); QLon=dat(:,1);QLat=dat(:,2);Qf(:,1)=dat(:,3);
dat=load('Inv1_Q_1_5'); Qf(:,2)=dat(:,3);
dat=load('Inv1_Q_3'); Qf(:,3)=dat(:,3);
dat=load('Inv1_Q_6');Qf(:,4)=dat(:,3);
dat=load('Inv1_Q_12');Qf(:,5)=dat(:,3);


Q0=0*QLat+NaN;
eta=0*QLat+NaN;
for i=1:numel(QLon)
    okay=find(~isnan(Qf(i,:)));
    if length(okay)>=4
        Qs=Qf(i,okay);
        fs=frequencies(okay);
        lQs=log(Qs)';
        lfs=log(fs);
        p=polyfit(lfs,lQs,1);
        eta(i)=p(1);
        Q0(i)=exp(p(2));
    end
end
make_contour(QLon,QLat,Q0,0.25);title('Q0')
make_contour(QLon,QLat,eta,0.25);title('eta')
dlmwrite('Q0', [QLon QLat Q0],'\t')
dlmwrite('eta', [QLon QLat eta],'\t')
