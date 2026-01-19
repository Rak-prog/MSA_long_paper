function out =  merge_networks(T, nets, left_right_merge)

    vn = string(T.Properties.VariableNames);

    if left_right_merge == false
        hemis = ["LH","RH"];
            
        out = table('Size',[height(T) numel(nets)*numel(hemis)], ...
            'VariableTypes', repmat("double",1,numel(nets)*numel(hemis)));
        colNames = strings(1, numel(nets)*numel(hemis));
        c = 0;
    
        for h = 1:numel(hemis)
            for k = 1:numel(nets)
                c = c + 1;
                colNames(c) = hemis(h) + "_" + nets(k);
        
                cols = contains(vn, "_" + hemis(h) + "_" + nets(k) + "_", 'IgnoreCase', false);
                out.(colNames(c)) = mean(T{:, cols}, 2, 'omitnan');
            end
        end
    elseif left_right_merge == true
        
        out = table('Size',[height(T) numel(nets)], ...
            'VariableTypes', repmat("double",1,numel(nets)), ...
            'VariableNames', cellstr(string(nets)));
        colNames = strings(1, numel(nets));
        c = 0;
        for k = 1:numel(nets)
            c = c + 1;
            colNames(c) = nets(k);
            % match either _LH_<net>_ or _RH_<net>_
            cols = contains(vn, "_LH_" + nets(k) + "_", 'IgnoreCase', false) | ...
                   contains(vn, "_RH_" + nets(k) + "_", 'IgnoreCase', false);

            out.(string(nets(k))) = mean(T{:, cols}, 2, 'omitnan');
        end

    end
    
    colsAllZero = all(out{:,:} == 0, 1);  
    out(:, colsAllZero) = [];
    out.Properties.VariableNames = cellstr(colNames);


end