%load('SimLab__1_05-Mar-2023.mat') % 1st dataset
load('SimLab_IE_2_19-Mar-2023.mat');
samplePoints = length(results(:,:,1,2));
testParameters = length(synFailProb);
sparsity_set1 = synFailProb;
numMetrics = 6;
meanMetrics = zeros(numMetrics,testParameters,samplePoints);
ieMetrics = zeros(numMetrics,testParameters,samplePoints);
for i = 1:samplePoints
    for j = 1:testParameters
    ieMetrics(:,j,i) = squeeze(mean(results{i,j,1,2}));
    end
end

%load('SimLab_IE_1_08-Mar-2023.mat'); % 1st dataset
for i = 1:samplePoints
        for j = 1:testParameters
    meanMetrics(:,j,i) = squeeze(mean(results{i,j,1,1}));
        end
end
% Metric indices are
% 1. False positive Rate. 2. False Negative Rate 3. True Correct Rate
% 4. Weighted Correct Rate 5. False Discovery Rate. 6. False Omission Rate
datasets = {ieMetrics, meanMetrics};
datanames = {'Accuracy','Recall','Precision','False Positive Rate','F1'};
figure('Renderer', 'painters', 'color','w','Position', [10 10 1600 900])
MetricIdx = [3 2 5 1 7];
ieMetrics(2,:,:) = 1 - ieMetrics(2,:,:); meanMetrics(2,:,:) = 1 - meanMetrics(2,:,:);
ieMetrics(5,:,:) = 1 - ieMetrics(5,:,:); meanMetrics(5,:,:) = 1 - meanMetrics(5,:,:);
ieMetrics(7,:,:) = 2*(ieMetrics(5,:,:).*ieMetrics(2,:,:))./(ieMetrics(5,:,:)+ieMetrics(2,:,:)); 
meanMetrics(7,:,:) = 2*(meanMetrics(5,:,:).*meanMetrics(2,:,:))./(meanMetrics(5,:,:)+meanMetrics(2,:,:)); 

data_indices = [2 3];
sp_idx = 0;
end_idx = 6;
for i = [2 3]%1:length(datanames)
    sp_idx = sp_idx + 1;
subplot(length(data_indices),3,2 + (sp_idx-1)*3)

%semilogx(N_obsCell./test_samples,[EOnlyMetrics(MetricIdx(i),:);BaseMetrics(MetricIdx(i),:)],'*-', 'LineWidth',2)
%set ( gca, 'xdir', 'reverse' )
%plot(test_samples,[EOnlyMetrics(MetricIdx(i),:);BaseMetrics(MetricIdx(i),:)],'*-', 'LineWidth',2)
plot(test_samples,squeeze(meanMetrics(MetricIdx(i),1:end_idx,:)),'*-', 'LineWidth',2)
hold on
%plot([0 50],[0 0], '--k', 'LineWidth',1)
colormap('jet')
%xlabel('Compression Ratio','FontSize',12,'FontWeight','bold')
ylim([0 1])
xlabel('Trials','FontSize',14,'FontWeight','bold')
ylabel(datanames{i},'FontSize',14,'FontWeight','bold')
%xticklabels(0:.2:1)
%legend({'E and I cell experiment','E cell only experiment'},'Location', 'southeast')
legend(num2str(sparsity_set1'),'Location', 'southeast')
title('Synaptic Failure  (I/E Constraints)','FontSize',14,'FontWeight','bold')

subplot(length(data_indices),3,1 + (sp_idx-1)*3)

%semilogx(N_obsCell./test_samples,100*((BaseMetrics(MetricIdx(i),:)) - (EOnlyMetrics(MetricIdx(i),:))),'*-', 'LineWidth',2)
%set ( gca, 'xdir', 'reverse' )
plot(test_samples,squeeze(ieMetrics(MetricIdx(i),1:end_idx,:)),'*-', 'LineWidth',2)
%plot(test_samples,squeeze(BaseMetrics(MetricIdx(i),:,:)),'*-', 'LineWidth',2)
hold on
%plot([0 50],[0 0], '--k', 'LineWidth',1)
colormap('jet')
ylim([0 1])
xlabel('Trials','FontSize',14,'FontWeight','bold')
ylabel([datanames{i}],'FontSize',14,'FontWeight','bold')
%xticklabels(0:.2:1)
%legend('E only - E and I','Location', 'southeast')
legend(num2str(sparsity_set1'),'Location', 'southeast')
%legend(num2str(round(20000./TotalCells(1:end-1)')),'Location', 'southeast')
title('Synaptic Failure (Baseline)','FontSize',14,'FontWeight','bold')

subplot(length(data_indices),3,3 + (sp_idx-1)*3)
diffMetricPlot = squeeze(meanMetrics(MetricIdx(i),1:end_idx,:))-squeeze(ieMetrics(MetricIdx(i),1:end_idx,:));
plot(test_samples,diffMetricPlot,'*-', 'LineWidth',2)
hold on
colormap('jet')
xlabel('Trials','FontSize',14,'FontWeight','bold')
ylabel([datanames{i}  ' Difference'],'FontSize',14,'FontWeight','bold')
legend(num2str(sparsity_set1'),'Location', 'northeast')
title('Synaptic Failure  (I/E - Baseline)','FontSize',14,'FontWeight','bold')


end