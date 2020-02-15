function [hfinal] = merge_posterior_figures_fromhandles(f_array,nplots,p)
    %%
    nfiles = length(f_array);
    for j = 1:nfiles
        % Your original figures
        h(j).fig = f_array(j);  
        figure_ax = findall(f_array(j),'type','axes');
        c = 1;
        for i = nplots:-1:1
            h(j).ax(i) = figure_ax(c);% subplot(1,nplots,i);
            c=c+1;
        end
    end
    if isempty(p)
        p1 = nfiles;
        p2 = nplots;
    else
        p1 = p(1);p2=p(2);
    end
    %%
    hfinal.fig = figure%('units','normalized','outerposition',[0 0 1 1]);
    hfinal.ax = gobjects(nplots*nfiles);
    [ha, pos] = tight_subplot(p1,p2);
    for ii = 1:nplots*nfiles
        hfinal.ax(ii) = ha(ii);%subplot(p1,p2,ii);
    end
    c = 1;
    for j = 1:nfiles
        for i =1:nplots
            hfinal.ax2(c,1) = copyobj(h(j).ax(i), hfinal.fig);
            c = c+1;
        end
    end
    for ii = 1:nplots*nfiles
        hfinal.ax2(ii).Position = hfinal.ax(ii).Position;
    end
        
    delete(hfinal.ax);
end