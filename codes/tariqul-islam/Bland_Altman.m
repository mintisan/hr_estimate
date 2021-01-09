function  []= Bland_Altman( Estimated_Values, Ground_Truth_Values, title__, x_label, y_label, xlim__, ylim__ )

    diff = Estimated_Values - Ground_Truth_Values;
    avg = (Estimated_Values+Ground_Truth_Values)./2;
    mn = mean(diff)
    sd = std(diff);
    val1 = mn - 1.96*sd;
    val2 = mn + 1.96*sd;
    
    LOA = [val1, val2]
    
    percentage = mean( diff>val1 & diff< val2 )
    
    figure;
    plot(avg,diff,'+');
    hold on;
    nbar = xlim__(1):xlim__(2);
    valbar = ones(1,xlim__(2)-xlim__(1)+1);
    plot(nbar,valbar*val1,'r--');
    text(xlim__(2)-2,val1-1.5,'\mu-1.96\sigma','HorizontalAlignment','right','Fontname','Times New Roman','FontSize', 14, 'Color', 'r');
    plot(nbar,valbar*val2,'r--');
    text(xlim__(2)-2,val2+2.2,'\mu+1.96\sigma','HorizontalAlignment','right','Fontname','Times New Roman','FontSize', 14, 'Color', 'r');
    plot(nbar,valbar*mn,'k-');
    text(xlim__(2)-2,mn+2.2,'\mu','HorizontalAlignment','right','Fontname','Times New Roman','FontSize', 14, 'Color', 'r');
    ylim(ylim__);
    xlim(xlim__);
    xlabel(x_label, 'Fontname','Times New Roman','FontSize', 16)
    ylabel(y_label, 'Fontname','Times New Roman','FontSize', 16)
    title(title__, 'fontweight','normal', 'Fontname','Times New Roman','FontSize', 17)
    

end
