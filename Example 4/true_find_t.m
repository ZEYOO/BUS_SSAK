
function t_true = true_find_t(p0,y_mcs)
Num = length(y_mcs);
y_sort = sort(y_mcs); 
t_index = Num*p0;
t_true = y_sort(t_index);
end