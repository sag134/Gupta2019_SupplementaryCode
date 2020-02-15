function  expt = get_expt_struct(mean_data,std_data)
    %take mean and standard deviation of data and create ptempest styled
    %expt struct for fitting.
    expt{1}.mean = mean_data;
    expt{1}.stdev = std_data;
    expt{1}.nsamples = ones(size(expt{1}.mean));
    expt{1}.weight = ones(size(expt{1}.mean));
    expt{1}.time = 1:size(expt{1}.mean,1);
end

