function [out_low, out_high] = get_outliers(data, subj_id, multiply)
data = table2array(data);
subj_id = table2array(subj_id);
low = mean(data) - multiply *  iqr(data,1);
high = mean(data) + multiply *  iqr(data,1);

mask_low = data < low;
[row_low, col_low] = find(mask_low);

mask_high = data > high;
[row_high, col_high] = find(mask_high);

out_low=subj_id(find(sum(mask_low,2)>5))

out_high=subj_id(find(sum(mask_high,2)>5))

end
