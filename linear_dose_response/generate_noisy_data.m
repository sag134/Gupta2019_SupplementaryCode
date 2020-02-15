addpath('../lib')
true_params = [1,-10,-10,-10,-2,-10,-10];
plt = 1;
noise_level = 0.1;
repeats = 10;
params = true_params;
example = 'linear_full_trajectory';
tmp_simulate_function = str2func(['simulate_',example]);
for i = 1:4
    simulate_function = @(params)tmp_simulate_function(i,[],[],params);
    [smooth_data,noisy_data,mean_data,std_data,noise,f]  = gen_noisy_data(params,simulate_function,noise_level,repeats,plt);    
    tmp = get_expt_struct(mean_data,std_data);
    expt{i} = tmp{1};
%save(['./',example,'.mat'],'expt')
%save([example,'_rawdata.mat'],'smooth_data','noisy_data','mean_data','std_data','noise','f');
end
   