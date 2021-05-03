
%%% this script handles the (first) inversion. It is written to behave like
%%% a function; inputs are re-named before use in a rather confusing way to
%%% preserve them for comparison or other uses. 

weights=ones(size(amp2));  %%% a placeholder in case weights are applicable
%%% Weighted inversions were validated on at least one occasion, but be
%%% sure to double check that weights work with the current implementation 
fprintf([ 'Before inversion 1. L2 norm=' num2str(mean(signal.^2)) '\n'])

dd=[signal.*weights; zeros(length(X),1)]; 
%%%%dd: [current log(amplitude) residuals
%%%      zeros attendant to smoothing matrix]
GG=[sonum2.*weights stanum2.*weights path_props.*weights; smooth_A];
%%% GG: [ (source ID---station ID---kernels) weighted
%%%       smoothing (penalty) function];

%%%% GG model1 ~=~ dd ... model1=GG^-1  dd 
model1=pinv(GG)*dd;
pred1=GG*model1;
misfit1=dd-pred1;
misfit1=misfit1(1:length(amp2));  %%% discard the smoothing zeros to retain only amplitude residuals
misfit1=misfit1./weights; %%% I'm not sure that I validated this step with non-constant "weights"
fprintf([ 'After inversion 1. L2 norm=' num2str(mean(misfit1(1:length(amp2)).^2)) '\n'])

dS=model1(1:nso2); %%% source terms relative to reference at this frequency
dR=model1(nso2+1:nso2+nsta2); %%% receiver terms
dA=model1(nso2+nsta2+1:end); %%% 1/Q deviations from constant 1/Q background

soterm=sonum2*dS; %%% does not include reference S term at this frequency
staterm=stanum2*dR;
A=dA+background_attenuation; %%% add average dataset-wide 1/Q back in

%%%% make sure that the forward calcs with S, R, and 1/Q models reproduce
%%%% residuals from the inversion.
signal_=signal-soterm-staterm-path_props*(A-background_attenuation); 
fprintf([' Checking reconstructed residual. L2 norm=' num2str(mean(signal_.^2)) '\n'])


%%% small 1/A (low Q) nodes are not reliable. What is reliable? 
%%% The median of all perturbations is robust (ought to be near 0, but it 
%%% doesn'e actually matter). Large (positive) A perturbations are more 
%%% robust than small A (as Q-->big).
%%% Take the range of positive A perturbations relative to median as the
%%% acceptable range of A on either side of median. 
%%% In other words, if the maximum A is 3*median(A), a cludge to set the
%%% minimum A to median(A)/3 is justified.
Arange=max(A)/median(A);
Amin=median(A)/Arange;
A(A<Amin)=Amin;
Q=1./A;

signal_=signal-soterm-staterm-path_props*(A-background_attenuation);
fprintf([ 'After Q limits. L2 norm=' num2str(mean(signal_.^2)) '\n'])


%%% with the low-Q nodes possibly set to the floor value Amin=median(A) /
%%% (max(A)/median(A)), run the inversion again to allow other nodes and
%%% the source and receiver terms to update to the limits.

dd=[signal_.*weights; zeros(length(X),1)]; 
%%%% [current log(amplitude) residuals
%%%   zeros attendant to smoothing matrix]

%%%%for a boundary function, make weights of very low-Q nodes much higher
%%%% while pulling back a bit on the rest of the smoothing. 
smooth_factor=0.5*(Q/median(Q)).^4;smooth_factor(smooth_factor<0.5)=0.5;


GG=[sonum2.*weights stanum2.*weights path_props.*weights; smooth_A.*smooth_factor];

model1A=pinv(GG)*dd;
ddS=model1A(1:nso2);
ddR=model1A(nso2+1:nso2+nsta2);
ddA=model1A(nso2+nsta2+1:end);

A=A+ddA;
Q=1./A;
dS=dS+ddS;
soterm=sonum2*dS;
dR=dR+ddR;
staterm=stanum2*dR;
make_contour(Lon,Lat,Q,0.25);colormap(flipud(jet));
title(['Inverted Q at ' num2str(f) ' Hz'])

signal_=signal-soterm-staterm-path_props*(A-background_attenuation);
fprintf([ 'After bounded-Q inversion. L2 norm=' num2str(mean(signal_.^2)) '\n'])













