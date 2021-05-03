function [c,h]=make_contour(X,Y,val,spacing,max_spacing)


        Xdim=max(X)-min(X);
        Ydim=max(Y)-min(Y);
        area=Xdim*Ydim;
        point_area=spacing^2;
        npts=area/point_area;
        if npts>1e6
                warning('Too many points! \n')

            point_area=area/5e5;
            newspacing=sqrt(point_area);
            warning([ num2str(spacing) ' m spacing too fine. Changing to ' num2str(newspacing) ' m '])
            disp([ num2str(spacing) ' m spacing too fine. Changing to ' num2str(newspacing) ' m '])
            spacing=newspacing;
        end
        

        if nargin==4;max_spacing=10*spacing;end
        XGrid=min(X)-spacing/2:spacing:max(X)+spacing/2;
        YGrid=min(Y)-spacing/2:spacing:max(Y)+spacing/2;
        [XGrid,YGrid]=meshgrid(XGrid,YGrid);
        ncontour=ceil(sqrt(length(X))*5)+10;
        F=scatteredInterpolant(X,Y,val,'natural','none');
        valGrid=F(XGrid,YGrid);
        
        
%         linspace( min(valGrid(~isnan(valGrid))), max(valGrid(~isnan(valGrid))), ncontour);
%         min(min(valGrid)),max(max(valGrid)),ncontour);

        for i=1:numel(valGrid)
             d=sqrt( (XGrid(i)-X).^2 + (YGrid(i)-Y).^2);
            if min(d)>max_spacing
%             weights=exp(- d.^2/(2*(2)).^2);
% % %             weights=exp(- d.^2/(2*(1.82*minAnomDepth+4.35)).^2);
%             weights=weights/sum(weights);
%             valGrid(i)=dot(weights,val);
%             else
% %             f=find( abs(XGrid(i)-X)<max_spacing & abs(YGrid(i)-Y)<max_spacing, 1 );
%             if isempty(f)
                valGrid(i)=NaN;
            end
            
%             if isnan(valGrid(i))
%                 pause
%             end
        end
        
         v=val(~isnan(val));
        mv=median(v);
        dv=v-mv;
        
        cmax=quantile(abs(dv),0.95);
        cmaxmax=quantile(abs(dv),0.999);lowcontour=mv-cmaxmax;highcontour=mv+cmaxmax;
        cmin=mv-cmax;
        cmax=mv+cmax;
        if cmin==cmax;cmin=min(val(~isnan(val)));cmax=max(val(~isnan(val)));end
                contours=linspace(lowcontour,highcontour,ncontour);
                
%                 min(val(~isnan(val))),max(val(~isnan(val))),ncontour);
        valGrid(valGrid>contours(end))=contours(end);
        valGrid(valGrid<contours(1))=contours(1);
        figure;
        [c,h]=contour(XGrid,YGrid,valGrid,	contours,'fill','on');
        
        colormap(jet);
        
        axis equal;axis tight
        xlim([min(X) max(X)])
        ylim([min(Y) max(Y)])
        hold on;
         ratio=(max(X)-min(X))/(max(Y)-min(Y));
        set(gcf,'position',[0 0 600*ratio*1.1 600])
        set(gca,'units','normalize','position',[0.04 0.03 0.94 0.95])
       colorbar
%         caxis([quantile(v,0.01) quantile(v,0.99)])
        caxis( [cmin cmax]);
%         colorbar%('Position',[0.1 0.5 0.025 0.4],'fontsize',12)
        
      
%      [~,~]=contour(XGrid,YGrid,valGrid,[ median( valGrid(~isnan(valGrid))) median( valGrid(~isnan(valGrid)))],  'fill','off' ,'color',[0.5 0.5 0.5]);