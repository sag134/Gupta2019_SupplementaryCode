function [mean_gf,std_gf] = get_gf_list(mu_range,b_range,file_prefix,file_suffix,n,step)
%%
    n_mu = length(mu_range); % different mu values tested
    n_b = length(b_range); % different b values tested 
    % we will calculate the mean and standard deviation of the likelihood
    % (goodness of fit) distribution for every combination of
    % hyperparameter
    mean_gf = zeros(n_mu,n_b);
    std_gf = zeros(n_mu,n_b);
    for i = 1:n_mu
        for j= 1:n_b
            i,j
            %load data for the particular hyperparameter set
            file_name = [file_prefix,'mu_',num2str(mu_range(i)),'_b_',num2str(b_range(j)),file_suffix];
            if exist(file_name, 'file') == 2
                data = load(file_name);
                ec = data.energy_chain;
                pc = data.params_chain;
                ec(1,1)
                % As long as the energy converged to a non infinite value
                if ec(1,1)<1e29
                    cfg = data.cfg;
                    N = size(ec,2);
                    c = 1;
                    gf = zeros(1,length(N-n:step:N));
                    for k = N-n:step:N
                        param = pc(1,:,k);
                        %Calculate the likelihood, i.e. the energy without
                        %the lasso penalty
                        gf(c) =energy_gaussian_wopenalty(param,cfg); 
                        c=c+1;
                    end
                    %calculate mean and standard deviation
                    mean_gf(i,j) = mean(gf);
                    std_gf(i,j) = std(gf);
                else
                    mean_gf(i,j) = NaN;
                    std_gf(i,j) = NaN;
                end
            else
                mean_gf(i,j) = NaN;
                std_gf(i,j) = NaN;
            end
        end
    end
end