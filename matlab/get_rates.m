function rates = get_rates(h, error, skip)
    h = h(skip+1:end);
    error = error(skip+1:end, :);
    
    rates = [log(h), ones(length(h),1)] \ log(error);
    rates = rates(1,:)';
end