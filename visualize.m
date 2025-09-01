function [h, mediansT0, mediansT1, iqrT0, iqrT1 , subj_outliers_tot] = visualize(dataT0, dataT1, titname, lim_y, labels, yname, subj_listT0, subj_listT1, BL_FUP, colT0, colT1)

h=figure('Color',[1 1 1]); hold on

nROIs = size(dataT0, 2);
nSubjects = size(dataT0, 1);

% Set spacing between ROI pairs
gap = 3;  % Try 3 or 4 for wider spacing
x0 = (0.5:gap:(nROIs-1)*gap+0.5);   % Positions for BL
x1 = x0 + 1;                        % Positions for FUP

% Subject lines
subj_outliers_tot = {};
for i = 1:nROIs 
    subj_out = {};

    for s = 1:nSubjects
        x = [x0(i), x1(i)];
        y = [dataT0(s,i), dataT1(s,i)];

        plot(x(1), y(1), [colT0,'o'], 'MarkerSize', 1, 'MarkerFaceColor', colT0);
        plot(x(2), y(2), [colT1,'o'], 'MarkerSize', 1, 'MarkerFaceColor', colT1);
        if BL_FUP == 1
            % --- CONDITION 1: Increase from BL to FUP ---
            if y(2) > y(1)
                plot(x, y, '-', 'Color', 'r', 'LineWidth', 0.8);
                subj_out{end+1} = subj_listT0{s}(1:6);
                if size(dataT0,2)<10
                    text(x(2)+0.2, y(2), subj_listT0{s}(1:6), 'FontSize', 4, 'Rotation', 0, 'Color', 'r');
                end
            else 
                plot(x, y, '-', 'Color', 'g', 'LineWidth', 0.8);
            end

            % --- CONDITION 2: Outlier (e.g., above +2 SD at FUP) ---
            % yFUP_col = dataT1(:,i);
            % meanFUP = mean(yFUP_col, 'omitnan');
            % if y(1) > meanBL + 2*stdBL || y(1) < meanBL - 2*stdBL
            %     text(x(1)+0.2, y(1), subj_listT1{s}, 'FontSize', 4, 'Rotation', 45, 'Color', 'b');
            % end
            % stdFUP = std(yFUP_col, 'omitnan');
            % 
            % if y(2) > meanFUP + 2*stdFUP || y(2) < meanFUP - 2*stdFUP
            %     text(x(2)+0.2, y(2), subj_list{s}, 'FontSize', 6, 'Rotation', 45, 'Color', 'r');
            % end

        elseif BL_FUP == 0

            % Extract current ROI column
            yBL_col = dataT0(:,i);
            yFUP_col = dataT1(:,i);

            % Compute mean and std for each timepoint USING MATLAB
            % DEFINITION OF OUTLIER ON BOXPLOT
            % meanBL = mean(yBL_col, 'omitnan');
            % stdBL = std(yBL_col, 'omitnan');
            % meanFUP = mean(yFUP_col, 'omitnan');
            % stdFUP = std(yFUP_col, 'omitnan');

            % --- Loop over BL subjects ---
            for s0 = 1:size(dataT0,1)
                y0 = dataT0(s0,i);
                if ~isnan(y0) && (y0 > prctile(yBL_col,75) + 1.5*iqr(yBL_col) || y0 < prctile(yBL_col,25) - 1.5*iqr(yBL_col) )
                    text(x0(i)+0.2, y0, subj_listT0{s0}, 'FontSize', 4, 'Rotation', 45, 'Color', 'k');
                end
            end

            % --- Loop over FUP subjects ---
            for s1 = 1:size(dataT1,1)
                y1 = dataT1(s1,i);
                if ~isnan(y1) && (y1 > prctile(yFUP_col,75) + 1.5*iqr(yFUP_col) || y1 < prctile(yFUP_col,25) - 1.5*iqr(yFUP_col))
                    text(x1(i)+0.2, y1, subj_listT1{s1}, 'FontSize', 4, 'Rotation', 45, 'Color', 'k');
                end
            end
        end
    end
    subj_outliers_tot{i} = subj_out;
end

% Boxplots
boxplot(dataT0, 'Positions', x0, 'Colors', colT0, 'Widths', 1);
boxplot(dataT1, 'Positions', x1, 'Colors', colT1, 'Widths', 1);
mediansT0 = mean(dataT0, 1); 
mediansT1 = mean(dataT1, 1);
iqrT0 = std(dataT0);
iqrT1 = std(dataT1);

if size(dataT0,2) == 3 & BL_FUP == 0 % do this only for the the volumes
 
    for i = 1:length(x0)
         text(x0(i), mediansT0(i), sprintf('%.2f', mediansT0(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 6,  'FontWeight', 'bold', 'Color', 'k');
         text(x1(i), mediansT1(i), sprintf('%.2f', mediansT1(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 6,  'FontWeight', 'bold', 'Color', 'k');
    end
end

% Axis formatting
midpoints = (x0 + x1) / 2;
set(gca, 'XTick', midpoints, 'XTickLabel', labels, 'XTickLabelRotation', 45);
xlim([min(x0)-1, max(x1)+1]);
ylabel(yname);
ylim(lim_y);
title(titname, 'Interpreter', 'none');

box off;

end