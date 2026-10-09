close all; clearvars; clc

%% import
%       Data already filtered with FastIca must be provided. If not use
%       also next section and adapt the code
place=1;
for i=30:40
    if i==30 || i==31 
        user=15;
    end
    if i==32 || i==33
        user=16;
    end
    if i==34 || i==35
        user=17;
    end
    if i==36 || i==37
        user=18;
    end
    if i==38 || i==39
        user=19;
    end
    if i==40
        user=20;
    end
    
    exp=int2str(i);
    user=int2str(user);
    
    file=strcat('ICA_exp',exp,'_user',user,'.mat');
    all(place).data_all=load(file);
    place=place+1;
end  
clear place exp file user i

%% ICA
    %use this if data provided have not been cleaned with fastICA yet.  
% for i =1:11
%     [fast,A,w]=fastica(all(i).data_all.data');
%     A(:,6)=0;
%     recovered_ICA=A*fast;
%     all(i).recovered_ICA=recovered_ICA';
%     clear A fast w recovered_ICA
% end
% clear i
    
%% Frequency

lowpass=load('firlow2.mat');

%apply the filter to the signal
for i=1:11
    all(i).recovered_ICA=all(i).data_all.recovered;
    all(i).filtered_recovered_ICA=filter(lowpass.Num,1,all(i).recovered_ICA); %
end

clear a b lowpass i

%% wavelet
dec_level=5;                    
mother_wav='db3';

%I divided in segments the signals on all 6 sensors
for i=1:11
    for j=1:length(all(i).data_all.label_data)
        start=all(i).data_all.label_data(j,4);
        stop=all(i).data_all.label_data(j,5);
        divided(i).segments(j).p_segment_sens1=all(i).filtered_recovered_ICA(start:stop,1);
        divided(i).segments(j).p_segment_sens2=all(i).filtered_recovered_ICA(start:stop,2);
        divided(i).segments(j).p_segment_sens3=all(i).filtered_recovered_ICA(start:stop,3);
        divided(i).segments(j).p_segment_sens4=all(i).filtered_recovered_ICA(start:stop,4);
        divided(i).segments(j).p_segment_sens5=all(i).filtered_recovered_ICA(start:stop,5);
        divided(i).segments(j).p_segment_sens6=all(i).filtered_recovered_ICA(start:stop,6);
    end
end
clear i j start stop

% perform a multi-level decomposition of the signal using Wavelet
for i=1:11
    for j=1:length(all(i).data_all.label_data)
        [C,L] = wavedec (divided(i).segments(j).p_segment_sens1, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens1.C=C;
        divided(i).wavelet_C_L(j).sens1.L=L;

        [C,L] = wavedec (divided(i).segments(j).p_segment_sens2, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens2.C=C;
        divided(i).wavelet_C_L(j).sens2.L=L;

        [C,L] = wavedec (divided(i).segments(j).p_segment_sens3, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens3.C=C;
        divided(i).wavelet_C_L(j).sens3.L=L;

        [C,L] = wavedec (divided(i).segments(j).p_segment_sens4, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens4.C=C;
        divided(i).wavelet_C_L(j).sens4.L=L;

        [C,L] = wavedec (divided(i).segments(j).p_segment_sens5, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens5.C=C;
        divided(i).wavelet_C_L(j).sens5.L=L;

        [C,L] = wavedec (divided(i).segments(j).p_segment_sens6, dec_level, mother_wav);
        divided(i).wavelet_C_L(j).sens6.C=C;
        divided(i).wavelet_C_L(j).sens6.L=L;
    end
end

for k=1:11
    for i=1:length(all(k).data_all.label_data)
        
        %reconstruct sens1
        rec_sens1=zeros(dec_level+1,length(divided(k).segments(i).p_segment_sens1));
            for j=1:dec_level
                rec_sens1(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens1.C,divided(k).wavelet_C_L(i).sens1.L,mother_wav,j);
            end
        rec_sens1(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens1.C,divided(k).wavelet_C_L(i).sens1.L,mother_wav,j);
        divided(k).Recons(i).sens1=rec_sens1;
        
        %reconstruct sens2
        rec_sens2=zeros(dec_level+1,length(divided(k).segments(i).p_segment_sens2));
            for j=1:dec_level
            rec_sens2(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens2.C,divided(k).wavelet_C_L(i).sens2.L,mother_wav,j);
            end
        rec_sens2(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens2.C,divided(k).wavelet_C_L(i).sens2.L,mother_wav,j);
        divided(k).Recons(i).sens2=rec_sens2;
        
        
        %reconstruct sens3
        rec_sens3=zeros(dec_level+2,length(divided(k).segments(i).p_segment_sens3));
            for j=1:dec_level
            rec_sens3(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens3.C,divided(k).wavelet_C_L(i).sens3.L,mother_wav,j);
            end
       rec_sens3(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens3.C,divided(k).wavelet_C_L(i).sens3.L,mother_wav,j);
        divided(k).Recons(i).sens3=rec_sens3;
        
        
        %reconstruct sens4
        rec_sens4=zeros(dec_level+2,length(divided(k).segments(i).p_segment_sens4));
            for j=1:dec_level
            rec_sens4(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens4.C,divided(k).wavelet_C_L(i).sens4.L,mother_wav,j);
            end
        rec_sens4(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens4.C,divided(k).wavelet_C_L(i).sens4.L,mother_wav,j);
        divided(k).Recons(i).sens4=rec_sens4;
        
        
        %reconstruct sens5
        rec_sens5=zeros(dec_level+2,length(divided(k).segments(i).p_segment_sens5));
            for j=1:dec_level
            rec_sens5(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens1.C,divided(k).wavelet_C_L(i).sens5.L,mother_wav,j);
            end
        rec_sens5(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens5.C,divided(k).wavelet_C_L(i).sens5.L,mother_wav,j);
        divided(k).Recons(i).sens5=rec_sens5;
        
        
        %reconstruct sens6
        rec_sens6=zeros(dec_level+2,length(divided(k).segments(i).p_segment_sens6));
            for j=1:dec_level
            rec_sens6(j,:)=wrcoef('d',divided(k).wavelet_C_L(i).sens6.C,divided(k).wavelet_C_L(i).sens6.L,mother_wav,j);
            end
        rec_sens6(dec_level+1,:)=wrcoef('a',divided(k).wavelet_C_L(i).sens6.C,divided(k).wavelet_C_L(i).sens6.L,mother_wav,j);
        divided(k).Recons(i).sens6=rec_sens6;
        
    end 
end
clear rec_sens1 rec_sens2 rec_sens3 rec_sens4 rec_sens5 rec_sens6 i j k L C start stop mother_wav 


%% feature extraction

%note: 
%          1:dec_level ->detail decomposition
%          dec_level+1   ->aproximation decomposition
levels=[1:dec_level+1];

for i=1:11
    for j=1:length(all(i).data_all.label_data)
        for k=levels
           
            %mean
            features(i).mean(j).sens1(k)=mean(divided(i).Recons(j).sens1(k,:));
            features(i).mean(j).sens2(k)=mean(divided(i).Recons(j).sens2(k,:));
            features(i).mean(j).sens3(k)=mean(divided(i).Recons(j).sens3(k,:));
            features(i).mean(j).sens4(k)=mean(divided(i).Recons(j).sens4(k,:));
            features(i).mean(j).sens5(k)=mean(divided(i).Recons(j).sens5(k,:));
            features(i).mean(j).sens6(k)=mean(divided(i).Recons(j).sens6(k,:));
            
            %Shannon Entropy
            features(i).E(j).sens1(k)=wentropy(divided(i).Recons(j).sens1(k,:),'shannon');
            features(i).E(j).sens2(k)=wentropy(divided(i).Recons(j).sens2(k,:),'shannon');
            features(i).E(j).sens3(k)=wentropy(divided(i).Recons(j).sens3(k,:),'shannon');
            features(i).E(j).sens4(k)=wentropy(divided(i).Recons(j).sens4(k,:),'shannon');
            features(i).E(j).sens5(k)=wentropy(divided(i).Recons(j).sens5(k,:),'shannon');
            features(i).E(j).sens6(k)=wentropy(divided(i).Recons(j).sens6(k,:),'shannon');
            
            %Standard Deviation
            features(i).std(j).sens1(k)=std(divided(i).Recons(j).sens1(k,:));
            features(i).std(j).sens2(k)=std(divided(i).Recons(j).sens2(k,:));
            features(i).std(j).sens3(k)=std(divided(i).Recons(j).sens3(k,:));
            features(i).std(j).sens4(k)=std(divided(i).Recons(j).sens4(k,:));
            features(i).std(j).sens5(k)=std(divided(i).Recons(j).sens5(k,:));
            features(i).std(j).sens6(k)=std(divided(i).Recons(j).sens6(k,:));
        end
    end
end


%save the activities in a vector
for k=1:11
    for i=1:length(all(k).data_all.label_data)
        features(k).activity(i)=double(all(k).data_all.label_data(i,3));
    end
end

clear i j k levels

%% create feature table

E=[]; S=[]; m=[];
Ent=[]; STD=[]; media=[]; 
for k=1:11
    for j=1:20
        E=[[features(k).E(j).sens1],[features(k).E(j).sens2],[features(k).E(j).sens3],[features(k).E(j).sens4],[features(k).E(j).sens5],[features(k).E(j).sens6]];
        Ent=[Ent;E];
        S=[[features(k).std(j).sens1],[features(k).std(j).sens2],[features(k).std(j).sens3],[features(k).std(j).sens4],[features(k).std(j).sens5],[features(k).std(j).sens6]];
        STD=[STD;S];
        m=[[features(k).mean(j).sens1],[features(k).mean(j).sens2],[features(k).mean(j).sens3],[features(k).mean(j).sens4],[features(k).mean(j).sens5],[features(k).mean(j).sens6]];
        media=[media;m];
    end
end
act=[];
activity=[];
for k=1:11
    act=[features(k).activity]';
    activity=[activity;act];
end

%normalize data
feature_table=zscore([STD,Ent,media]);
%remove outliers 
feature_table(abs(feature_table)>3)=0;

clear act E S m p k j Ent media S STD
%% Kruskal

data.X=feature_table';
data.y=activity';
data.dim=size(data.X,1);
data.num_data=size(data.X,2);

col_names=num2cell(1:size(feature_table,2));

% get Kruskal-Wallis scores

for i = 1:size(data.X,1)
    [P, ANOVA] = kruskalwallis(data.X(i,:),data.y,'off');
    
    P_values(i) = P; 
    h_values(i) = ANOVA{2,5};
    
    rank{i,1} = col_names(i);
    rank{i,2} = h_values(i);
end

[sorted_ranks,sorted_idx] = sort([rank{:,2}],2,'descend');

display = [];
for i=1:size(data.X,1)
   display = [display sprintf('%s ---> %.2f\n', int2str(rank{sorted_idx(i),1}{1}), rank{sorted_idx(i),2})];
end

disp(display) 

f=[];
for i=1:69
    f=[f feature_table(:,sorted_idx(i))];
end

feature_table=f;

clear f ANOVA col_names data display h_values i P P_values rank sorted_idx sorted_ranks

%% PCA
[coeff,Projected_data_PCA,~,~,explained,~] = pca(feature_table);
Projected_data_PCA(:,1)=0.5.*Projected_data_PCA(:,1);
figure
bar(explained(1:10))
title('Explained variance PCA')

%MDS
D = pdist(feature_table,'euclidean'); %other interesting: Cityblock, Cosine, Correlation

%non classical
Projected_data_MDS = mdscale(D,2,'criterion','strain');
Projected_data_MDS(:,1)=0.5.*Projected_data_MDS(:,1);

clear D coeff explained

%% define macrogroups

main3act=zeros(length(activity),1);
for i=1:length(activity)
    if activity(i)<=3
        main3act(i)=1;
    else if (activity(i)>3 && activity(i)<=6)
            main3act(i)=2;
        else 
            main3act(i)=3;
        end
    end
end
clear i
Projected_data_MDS_activity=[Projected_data_MDS,activity];

% create a table with only the dynamic activities
Projected_data_MDS_dynamics=Projected_data_MDS_activity(Projected_data_MDS_activity(:,end)==1|Projected_data_MDS_activity(:,end)==2|Projected_data_MDS_activity(:,end)==3,:);

%class vector for dynamic activities
dynamic_act=Projected_data_MDS_dynamics(:,end);
Projected_data_MDS_dynamics=Projected_data_MDS_dynamics(:,1:2);


%% plot all activities

figure
subplot (1,2,1)
scatter(Projected_data_PCA(:,1),Projected_data_PCA(:,2),20,activity)
title('PCA')
subplot(1,2,2)
scatter(Projected_data_MDS(:,1),Projected_data_MDS(:,2),30,activity)
title('MDS')

%% plot main3activities (static, dynamic, transition)
figure
for i=1:length(main3act)
    if main3act(i)==1
        plot(Projected_data_MDS(i,1),Projected_data_MDS(i,2),'Marker','.','MarkerSize',15, 'MarkerEdgeColor','b')
    else if main3act(i)==2
            plot(Projected_data_MDS(i,1),Projected_data_MDS(i,2),'Marker','.','MarkerSize',15,'MarkerEdgeColor','g')
        else if main3act(i)==3
            plot(Projected_data_MDS(i,1),Projected_data_MDS(i,2),'Marker','.','MarkerSize',15,'MarkerEdgeColor','r')
            end
        end
    end
    hold on
end
title('Static, Dynamic and Transition')

%% plot only dynamic activities
%plot WALKING-blue, WALKING UPSTAIRS-green, WALKING DOWNSTAIRS-red
figure
for i=1:length(activity)
    if Projected_data_MDS_activity(i,end)==1
        plot(Projected_data_MDS_activity(i,1),Projected_data_MDS_activity(i,2),'Marker','.','MarkerSize',15,'MarkerEdgeColor','b')
    else if Projected_data_MDS_activity(i,end)==2
            plot(Projected_data_MDS_activity(i,1),Projected_data_MDS_activity(i,2),'Marker','.','MarkerSize',15,'MarkerEdgeColor','g')
        else if Projected_data_MDS_activity(i,end)==3
            plot(Projected_data_MDS_activity(i,1),Projected_data_MDS_activity(i,2),'Marker','.','MarkerSize',15,'MarkerEdgeColor','r')
            end
        end
    end
    hold on
end
title('Dynamic Activities')
%% ------------------------------------------------------------------------
   %CLUSTERING

% CLUSTERING  k-means (k medoids is exactly the same kmeans->kmedoids)

%Static, Dynamic and Transition
[idx_main_K,C_main_K] = kmeans(Projected_data_MDS(:,1:2),3,'Replicates',5);
figure
plot(Projected_data_MDS(idx_main_K==1,1),Projected_data_MDS(idx_main_K==1,2),'r.','MarkerSize',15)
hold on
plot(Projected_data_MDS(idx_main_K==2,1),Projected_data_MDS(idx_main_K==2,2),'b.','MarkerSize',15)
hold on
plot(Projected_data_MDS(idx_main_K==3,1),Projected_data_MDS(idx_main_K==3,2),'g.','MarkerSize',15)
plot(C_main_K(:,1),C_main_K(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3')
hold off
title('CLUSTERING: Static, Dynamic and Transition')

% dynamic activities
[idx_dyn_K,C_dyn_K] = kmeans(Projected_data_MDS_dynamics(:,1:2),3);
figure;
plot(Projected_data_MDS_dynamics(idx_dyn_K==1,1),Projected_data_MDS_dynamics(idx_dyn_K==1,2),'r.','MarkerSize',10)
hold on
plot(Projected_data_MDS_dynamics(idx_dyn_K==2,1),Projected_data_MDS_dynamics(idx_dyn_K==2,2),'g.','MarkerSize',10)
hold on
plot(Projected_data_MDS_dynamics(idx_dyn_K==3,1),Projected_data_MDS_dynamics(idx_dyn_K==3,2),'b.','MarkerSize',10)
plot(C_dyn_K(:,1),C_dyn_K(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3')
hold off
title('K-means clustering')

%ALL THE ACTIVITIES
figure
[idx_all_K,C_all_K] = kmeans(Projected_data_MDS(:,[1,2]),12,'Replicates',5);
plot(Projected_data_MDS(idx_all_K==1,1),Projected_data_MDS(idx_all_K==1,2),'.','MarkerSize',10,'color',[0 0 0])
hold on
plot(Projected_data_MDS(idx_all_K==2,1),Projected_data_MDS(idx_all_K==2,2),'.','MarkerSize',10,'color',[0 0 0.5])
hold on
plot(Projected_data_MDS(idx_all_K==3,1),Projected_data_MDS(idx_all_K==3,2),'.','MarkerSize',10,'color',[0 0.5 0])
hold on
plot(Projected_data_MDS(idx_all_K==4,1),Projected_data_MDS(idx_all_K==4,2),'.','MarkerSize',10,'color',[0.5 0 0])
hold on
plot(Projected_data_MDS(idx_all_K==5,1),Projected_data_MDS(idx_all_K==5,2),'.','MarkerSize',10,'color',[0 0 1])
hold on
plot(Projected_data_MDS(idx_all_K==6,1),Projected_data_MDS(idx_all_K==6,2),'.','MarkerSize',10,'color',[0 1 0])
hold on
plot(Projected_data_MDS(idx_all_K==7,1),Projected_data_MDS(idx_all_K==7,2),'.','MarkerSize',10,'color',[1 0 0])
hold on
plot(Projected_data_MDS(idx_all_K==8,1),Projected_data_MDS(idx_all_K==8,2),'.','MarkerSize',10,'color',[0 0.5 1])
hold on
plot(Projected_data_MDS(idx_all_K==9,1),Projected_data_MDS(idx_all_K==9,2),'.','MarkerSize',10,'color',[0 1 0.5])
hold on
plot(Projected_data_MDS(idx_all_K==10,1),Projected_data_MDS(idx_all_K==10,2),'.','MarkerSize',10,'color',[0.5 1 0])
hold on
plot(Projected_data_MDS(idx_all_K==11,1),Projected_data_MDS(idx_all_K==11,2),'.','MarkerSize',10,'color',[1 1 0])
hold on
plot(Projected_data_MDS(idx_all_K==12,1),Projected_data_MDS(idx_all_K==12,2),'.','MarkerSize',10,'color',[0 1 1])
hold on
plot(C_all_K(:,1),C_all_K(:,2),'kx','MarkerSize',7,'LineWidth',3)
title('K-means clustering')

%% k medoids
[idx_main_Kmed,C_main_Kmed] = kmedoids(Projected_data_MDS(:,1:2),3,'Replicates',5);
[idx_dyn_Kmed,C_dyn_Kmed] = kmedoids(Projected_data_MDS_dynamics(:,1:2),3);
[idx_all_Kmed,C_all_Kmed] = kmedoids(Projected_data_MDS(:,[1,2]),12,'Replicates',5);

%% hyerarchical clustering

%3main
T_main=clusterdata(Projected_data_MDS,'Maxclust',3);
figure
scatter(Projected_data_MDS(:,1),Projected_data_MDS(:,2),30,T_main)

%dynamic
T_dyn=clusterdata(Projected_data_MDS_dynamics,'distance','euclidean','linkage','ward','Maxclust',3);
figure
scatter(Projected_data_MDS_dynamics(:,1),Projected_data_MDS_dynamics(:,2),20,T_dyn)

%all
T_all=clusterdata(Projected_data_MDS,'distance','euclidean','linkage','ward','Maxclust',12);
figure
scatter(Projected_data_MDS(:,1),Projected_data_MDS(:,2),20,T_all)
title('Hierarchical clustering')

%% internal evaluation clustering: silhouette coefficient
  
myfunc1 = @(X,K)(kmedoids(Projected_data_MDS(:,1:2),K,'Replicates',5));
myfunc2 = @(X,K)(kmedoids(Projected_data_MDS_dynamics(:,1:2),K,'Replicates',5));

%3 MAIN
    %k means
    Kmeans_eva = evalclusters(Projected_data_MDS,'kmeans','silhouette','klist',[1:12]);
    Silh_main_kmeans=Kmeans_eva.CriterionValues(3); 
    %kmedoids
    Kmedoids_eva = evalclusters(Projected_data_MDS,myfunc1,'silhouette','klist',[1:12]);
    Silh_main_medoids=Kmedoids_eva.CriterionValues(3);
    %hierarchical clustering
    Hier_eva=evalclusters(Projected_data_MDS,'linkage','silhouette','klist',[1:12]);
    Silh_main_hier=Hier_eva.CriterionValues(3); 
   
%DYNAMICS
    %k means
    Kmeans_eva_dyn = evalclusters(Projected_data_MDS_dynamics,'kmeans','silhouette','klist',[1:6]);
    Silh_dyn_kmeans=Kmeans_eva_dyn.CriterionValues(3);    
    %kmedoids
    Kmedoids_eva_dyn = evalclusters(Projected_data_MDS_dynamics,myfunc2,'silhouette','klist',[1:6]);
    Silh_dyn_medoids=Kmedoids_eva_dyn.CriterionValues(3);
    %hierarchical
    Hier_eva_dyn=evalclusters(Projected_data_MDS_dynamics,'linkage','silhouette','klist',[1:12]);
    Silh_dyn_hier=Hier_eva_dyn.CriterionValues(3);
    
%ALL
    %k means
    Silh_all_kmeans=Kmeans_eva.CriterionValues(12); 
    %kmedoids
    Silh_all_medoids=Kmedoids_eva.CriterionValues(12);
    %hierarchical
    Silh_all_hier=Hier_eva.CriterionValues(12); 


%% silhouette plots

figure
subplot(3,3,1)
silhouette(Projected_data_MDS,idx_main_K)
title('K-means main activities')

subplot(3,3,2)
silhouette(Projected_data_MDS_dynamics,idx_dyn_K)
title('K-means dynamic activities')

subplot(3,3,3)
silhouette(Projected_data_MDS,idx_all_K)
title('K-means all activities')

subplot(3,3,4)
silhouette(Projected_data_MDS,idx_main_Kmed)
title('K-medoids main activities')

subplot(3,3,5)
silhouette(Projected_data_MDS_dynamics,idx_dyn_Kmed)
title('K-medoids dynamic activities')

subplot(3,3,6)
silhouette(Projected_data_MDS,idx_all_Kmed)
title('K-medoids all activities')

subplot(3,3,7)
silhouette(Projected_data_MDS,T_main)
title('Hierarchical main activities')

subplot(3,3,8)
silhouette(Projected_data_MDS_dynamics,T_dyn)
title('Hierarchical dynamic activities')

subplot(3,3,9)
silhouette(Projected_data_MDS,T_all)
title('Hierarchical all activities')

%% external evaluation-------------------------------------------------------

%3 main
 Acc_main_Kmean=AccMeasure(main3act,idx_main_K);
 Acc_main_Kmed=AccMeasure(main3act,idx_main_Kmed);
 Acc_main_H=AccMeasure(main3act,T_main);

%dynamics
 Acc_dyn_Kmean=AccMeasure(dynamic_act,idx_dyn_K);
 Acc_dyn_kmed=AccMeasure(dynamic_act,idx_dyn_Kmed)
 Acc_dyn_H=AccMeasure(dynamic_act,T_dyn);

%all
  %I do it because AccMeasure doesn't work with 12 classes
 T_all(T_all==12)=11;
 activity(activity==12)=11;
 idx_all_K(idx_all_K==12)=11;
 idx_all_Kmed(idx_all_Kmed==12)=11;  
    
Acc_all_Kmean=AccMeasure(activity, idx_all_K);
Acc_all_Kmed=AccMeasure(activity, idx_all_Kmed);
Acc_all_H=AccMeasure(activity,T_all);




