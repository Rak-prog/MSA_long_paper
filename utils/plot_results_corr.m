function h = plot_results_corr(T, titname)
    
    mat = table2array(T(1:end,2:end));
    labs_name = T.Properties.VariableNames(3:2:end);
    labs_y =  strrep(T.Labels, '_', '');
        
    idx = find(contains(labs_y,'Et')); 
    if ~isempty(idx); labs_y(idx) = {'Age'}; end
    idx = find(contains(labs_y,'Sesso')); 
    if ~isempty(idx); labs_y(idx) = {'Sex'}; end
    idx = find(contains(labs_y,'SottotipoMotorio')); 
    if ~isempty(idx); labs_y(idx) = {'Phenotype'}; end
    idx = find(contains(labs_y,'MCI_s__1_no_0_')); 
    if ~isempty(idx); labs_y(idx) = {'MCI'}; end
    idx = find(contains(labs_y,'Time_T1')); 
    if ~isempty(idx); labs_y(idx) = {'Time'}; end
    idx = find(contains(labs_y,'ROI:TimeT1')); 
    if ~isempty(idx); labs_y(idx) = {'ROIxTime'}; end
    

    labels_colored = cell(numel(labs_name),1);

    if size(T,2) == 63
        for ii=1:numel(labs_name)
            if ii <= 10
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', [95, 95, 95]/255, labs_name{ii});
            elseif ii >= 11 && ii <= 15
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', [55, 126, 184]/255, labs_name{ii});
            elseif ii >= 16 && ii <= 22
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', [0, 255, 255]/255, labs_name{ii});
            elseif ii >= 23 && ii <= 26
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', [152, 78, 163]/255, labs_name{ii});
            elseif ii >= 27 && ii <= 28
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', 	[166, 86, 40]/255, labs_name{ii});
            elseif ii >= 29 && ii <= 30
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', 	[247, 129, 191]/255, labs_name{ii});
            elseif ii >= 31 
                labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s',[153, 153, 153]/255, labs_name{ii});
            %else
            %    labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', [255, 127, 0]/255, labs_name{ii});
            end
        end
    elseif size(T,2) == 15
        cols = [[95, 95, 95]; [55, 126, 184];[0, 255, 255];[152, 78, 163];[166, 86, 40];	[247, 129, 191];[153, 153, 153]];
        for ii=1:numel(labs_name)
            labels_colored{ii} = sprintf('\\color[rgb]{%f, %f, %f}%s', cols(ii,:)/255, labs_name{ii});
        end
    end

    p_value = mat(:,2:2:end);
    estimates = mat(:,1:2:end)
    if p_value(:) < 0
        error('p values are not correct')
    end
    
        % Define thresholds
    th1 = 0.001;
    th2 = 0.01;
    th3 = 0.05;
    
    % Initialize a matrix of discrete categories:
    % 1 = green (< th1), 2 = orange (between th1 and th2), 3 = red (> th2)
    M_cat = zeros(size(p_value));
    M_cat(p_value < th1) = 1;
    M_cat(p_value >= th1 & p_value <= th2) = 2;
    M_cat(p_value > th2 & p_value <= th3) = 3;
    M_cat(p_value > th3) = 4;
    
    % Define colormap for these 3 categories: green, orange, red
    cmap = [0 1 0;    % green
            0.3 0.68 0.61 % dark green
            1 1 0; %  yellow
            1 0 0];   % red
    
    % Plot using imagesc
    h=figure('Color',[1 1 1])
    imagesc(M_cat);
    colormap(cmap); hold on

    % Fix color axis to have 3 discrete colors exactly
    caxis([0.5 4.5]);
    
    [nRows, nCols] = size(M_cat);

    % Draw vertical lines
    for c = 0.5 : 1 : nCols+0.5
        plot([c c], [0.5 nRows+0.5], 'k-', 'LineWidth', 1);
    end

    % Draw horizontal lines
    for r = 0.5 : 1 : nRows+0.5
        plot([0.5 nCols+0.5], [r r], 'k-', 'LineWidth', 1);
    end

    colorbar('Ticks',1:4,'TickLabels',{'<0.001','0.001-0.01','0.01-0.05','>0.05'});
    set(gca,'XTick',1:1:size(p_value,2),'XTickLabel',labels_colored,'FontSize',8);
    set(gca,'YTick',1:1:size(p_value,1),'YTickLabel',labs_y,'FontSize',8);
    xtickangle(45)
    grid on

    title(titname); box off
    
    numel(find((M_cat==1)))
    numel(find((M_cat==2)))
    numel(find((M_cat==3)))
    numel(find((M_cat==4)))
end