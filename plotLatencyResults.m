 %latency sd results
load('SimLabLatsd_5_01-Aug-2022.mat')
samplePoints = length(results(:,:,1,2));
testParameters = length(latency_sd);
latency_sd_set1 = latency_sd;
numMetrics = 6;
meanMetrics = zeros(numMetrics,testParameters,samplePoints);
sdMetrics = zeros(numMetrics,testParameters,samplePoints);
for i = 1:samplePoints
    for j = 1:testParameters
    sdMetrics(:,j,i) = squeeze(mean(results{i,j,1,1}));
    end
end

 %latency mean results
load('SimLab_latmean_1_01-Aug-2022.mat');
for i = 1:samplePoints
        for j = 1:testParameters
    meanMetrics(:,j,i) = squeeze(mean(results{i,j,1,1}));
        end
end
% Metric indices are
% 1. False positive Rate. 2. False Negative Rate 3. True Correct Rate
% 4. Weighted Correct Rate 5. False Discovery Rate. 6. False Omission Rate
datasets = {sdMetrics, meanMetrics};
datanames = {'Accuracy','Recall','Precision','False Positive Rate','F1'};
figure('Renderer', 'painters', 'color','w','Position', [10 10 1600 900])
MetricIdx = [3 2 5 1 7];
sdMetrics(2,:,:) = 1 - sdMetrics(2,:,:); meanMetrics(2,:,:) = 1 - meanMetrics(2,:,:);
sdMetrics(5,:,:) = 1 - sdMetrics(5,:,:); meanMetrics(5,:,:) = 1 - meanMetrics(5,:,:);
sdMetrics(7,:,:) = 2*(sdMetrics(5,:,:).*sdMetrics(2,:,:))./(sdMetrics(5,:,:)+sdMetrics(2,:,:)); 
meanMetrics(7,:,:) = 2*(meanMetrics(5,:,:).*meanMetrics(2,:,:))./(meanMetrics(5,:,:)+meanMetrics(2,:,:)); 

data_indices = [2 3];
sp_idx = 0;
for i = [2 3]%1:length(datanames)
    sp_idx = sp_idx + 1;
subplot(2,length(data_indices),sp_idx)

%semilogx(N_obsCell./test_samples,[EOnlyMetrics(MetricIdx(i),:);BaseMetrics(MetricIdx(i),:)],'*-', 'LineWidth',2)
%set ( gca, 'xdir', 'reverse' )
%plot(test_samples,[EOnlyMetrics(MetricIdx(i),:);BaseMetrics(MetricIdx(i),:)],'*-', 'LineWidth',2)
plot(test_samples,squeeze(meanMetrics(MetricIdx(i),:,:)),'*-', 'LineWidth',2)
hold on
%plot([0 50],[0 0], '--k', 'LineWidth',1)
colormap('default')
%xlabel('Compression Ratio','FontSize',12,'FontWeight','bold')
xlabel('Trials','FontSize',14,'FontWeight','bold')
ylabel(datanames{i},'FontSize',14,'FontWeight','bold')
%xticklabels(0:.2:1)
%legend({'E and I cell experiment','E cell only experiment'},'Location', 'southeast')
legend(num2str(latency_mean'),'Location', 'southeast')
title('Latency mean (ms)','FontSize',14,'FontWeight','bold')

subplot(2,length(data_indices),sp_idx+length(data_indices))

%semilogx(N_obsCell./test_samples,100*((BaseMetrics(MetricIdx(i),:)) - (EOnlyMetrics(MetricIdx(i),:))),'*-', 'LineWidth',2)
%set ( gca, 'xdir', 'reverse' )
plot(test_samples,squeeze(sdMetrics(MetricIdx(i),:,:)),'*-', 'LineWidth',2)
%plot(test_samples,squeeze(BaseMetrics(MetricIdx(i),:,:)),'*-', 'LineWidth',2)
hold on
%plot([0 50],[0 0], '--k', 'LineWidth',1)
colormap('default')
xlabel('Trials','FontSize',14,'FontWeight','bold')
ylabel([datanames{i}],'FontSize',14,'FontWeight','bold')
%xticklabels(0:.2:1)
%legend('E only - E and I','Location', 'southeast')
legend(num2str(latency_sd_set1'),'Location', 'southeast')
%legend(num2str(round(20000./TotalCells(1:end-1)')),'Location', 'southeast')
title('Latency Standard Deviation (ms)','FontSize',14,'FontWeight','bold')

end