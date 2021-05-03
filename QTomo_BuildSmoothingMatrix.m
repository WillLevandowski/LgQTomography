 %%% 
    smoothing=zeros(length(X));
    smoothing_factor=smooth_coefficient*spacing^2/100^2; 
    %%% "smooth_coefficient" is defined in QTomo_SetParams
    %%% It and the exponent 0.8, below, were chosen by Tikhonov 
    %%% regularization of each frequency's inversion. These
        all_dist=zeros(length(X));
        for i=1:length(X);  all_dist(:,i)=sqrt((X-X(i)).^2+(Y-Y(i)).^2);      end
    
    %%% this smoothing function works well, subjectively. A Gaussian would
    %%% probably be fine, too. 
        for i=1:length(X)
            distances=all_dist(:,i);%sqrt((X-X(i)).^2+(Y-Y(i)).^2);
            smoothing(i,:)=smoothing_factor*f^smooth_f_exponent * (distances.^-0.5);%.* (0.5.^(d'/1000)) ;%/ sqrt(hits(i)+1);
            smoothing(i,i)=0;
            smoothing(i,i)=-sum(smoothing(i,:)); 
        end
    smooth_A=[zeros(length(X),length(stations2)+length(sources2)) smoothing];
 