addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/');
addpath('/shared2/LabUserFiles/Sanjana_Gupta/Original/ptempest/core/distr/');
addpath('../lib');
%%
fs = 12;
% p25 is activation penalty and p26 is IkB penalty
clr_fr = [0.8,0.8,1];
tcounter=0;
repeat_index = 1;%We are plotting covariation from just one representative repeat
for trajectory_number = [2] % and one representative trajectory
    tcounter=tcounter+1;
    LB = 5e4;
    UB = 9.9e5;
    num_repeats = length(repeat_index);
    dim = 26;
    num_pts = UB-LB+1;
    par_matrix = zeros(num_pts,dim,num_repeats);
    tmp = cell(1,num_repeats);
    % get data
    mu = -25;
    b = 2;
    for repeat = 1:num_repeats
        display(['trajectory_number: ',num2str(trajectory_number),' repeat: ',num2str(repeat_index(repeat))]);
        data = load(['../SingleCellNFkB_pulse/SingleCellNFkB_continue_chain_reduced_model_fixedstart_withlasso_mu_',num2str(mu),'_b_',num2str(b),'_trajectory_',num2str(trajectory_number),'_repeat_',num2str(repeat_index(repeat)),'_progress990000.mat']);    
        params = data.params_chain(1,:,LB:UB);
        tmp{repeat} = reshape(params,dim,num_pts);
        par_matrix(:,:,repeat) = tmp{repeat}';
    end
end
%%
f2 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);

p1 = par_matrix(:,13,:)+par_matrix(:,26,:); % IkB module
p2 = par_matrix(:,14,:)+par_matrix(:,26,:); % IkB module

s = binscatter(p1(:),p2(:),100);

xlim([-20 10]);ylim([-20 10])

hold on
plot([-20 10],[-20,10],'--k','LineWidth',2)
xlabel('Exit rate of free nuclear IkB')
ylabel('Exit rate of nuclear IkB-NFkB complex');


subplot(2,2,2);

p1 = par_matrix(:,11,:)+par_matrix(:,25,:); %Activation module
p2 = par_matrix(:,14,:)+par_matrix(:,26,:); % IkB module


binscatter(p1(:),p2(:),100)


xlim([-20 10]);ylim([-20 10])

hold on
plot([-20 10],[-20,10],'--k','LineWidth',2)
xlabel('Exit rate of free nuclear NFkB')
ylabel('Exit rate of nuclear IkB-NFkB complex');

subplot(2,2,3);

p1 = par_matrix(:,18,:)+par_matrix(:,25,:); %Activation module
p2 = par_matrix(:,16,:)+par_matrix(:,26,:); % IkB module


binscatter(p1(:),p2(:),100)

xlim([-20 10]);ylim([-20 10])

hold on
plot([-20 10],[-20,10],'--k','LineWidth',2)
xlabel('IKK mediated degradation of bound IkB')
ylabel('Basal degradation of bound IkB');

subplot(2,2,4);

p1 = par_matrix(:,17,:)+par_matrix(:,25,:); %Activation module
p2 = par_matrix(:,15,:)+par_matrix(:,26,:); % IkB module


binscatter(p1(:),p2(:),100)

xlim([-20 10]);ylim([-20 10])

hold on
plot([-20 10],[-20,10],'--k','LineWidth',2)
xlabel('IKK mediated degradation of free IkB')
ylabel('Basal degradation of free IkB');
%%
%toPDF(f2,'pdf_figures/FigureS5.pdf');