filename='ux_04_sh.txt';
y=dlmread(filename);
n=length(y);
wind=500;
tstep=5;
thr_fr=0.5;
crp_emb_dimm=1;
crp_delay=1;
crp_epsilon=0.4; %0.1 - Marwan default
logfile='ux_04_sh_vg_500_5.txt';
graph_type=2;%1 - crp, 2 - visibility
%graph_type=2;%visibility
% 20130601 Graph Entropy
Q=5;
types=round((1:wind+1)/(wind/(Q-1))+1);% вектор классов узлов (каждый из узлов должен быть отнесен к какому-то классу)
Nrunmax=3;% кол-во запусков процедуры
Nd=10;% кількість градацій відстані для ентропії відстаней
   % l=logspace(log(d_min)/log(10),log(Dmax)/log(10),Nd);
   %l=linspace(-14,55,30);
 % l=linspace(0.2,9,Nd);
skip_complex_measures=1;
% end 20130601 Graph Entropy
mas_diam=[];
mas_ld=[];
mas_sl=[];
mas_ncc=[];
mas_conn=[];
mas_simplgr=[];
mas_reggr=[];
mas_complgr=[];
mas_eulgr=[];
mas_treegr=[];
mas_graphic_gr=[];
mas_degr=[];
mas_maxdegr=[];
mas_mindegr=[];
mas_meandegr=[];
mas_avrndegr=[];
mas_maxavrndegr=[];
mas_minavrndegr=[];
mas_meanavrndegr=[];
% massivy:
% mas_max_closeness,mas_min_node_betwenness,mas_max_edge_betwenness,mas_max_weighted_clust_coeff,mas_max_vertex_eccentricity
mas_clust=[];
mas_average_degree=[];
mas_max_closeness=[];
mas_min_closeness=[];
mas_mean_closeness=[];
mas_min_node_betwenness=[];
mas_max_node_betwenness=[];
mas_mean_node_betwenness=[];
mas_min_edge_betwenness=[];
mas_max_edge_betwenness=[];
mas_mean_edge_betwenness=[];
mas_clust2=[];
mas_clust3=[];
    %[c1,c2,c3]=clust_coeff()
mas_min_weighted_clust_coeff=[];
mas_max_weighted_clust_coeff=[];
mas_mean_weighted_clust_coeff=[];
mas_ave_path_length=[];
mas_smooth_diameter=[];
mas_closeness=[];
mas_max_vertex_eccentricity=[];% max distance to any other node
mas_min_vertex_eccentricity=[];
mas_mean_vertex_eccentricity=[];
mas_theta=[];
mas_EntropyQ=[];
mas_Sigma_0=[];
mas_Sigma_d=[];
    % графика
    %draw_circ_graph
    %dot_matrix_plot()
    %radial_plot()
    % 20130703 убираются меры: self loops,num conn comp,connectivity,simple
    % graph\tregular graph\tcomplete graph\teulerian
    % graph\ttree\tgraphicбmin node degree\tmean node degree, min avg node
    % degree\tmean avg node degree,avg_degr\tmin(clsnss),\tmean(clsnss), mean(nb)
    %\\tmin(min((eb))\tmean(mean((eb))\tc2\\tmin(wcc)\tmean(wcc)\tapl\tsm_diam
    %\tmin(vert_ecc)\tmean(vert_ecc)\tEntropyQ\tTheta\tSigma_0\tSigma_d
fp=fopen(logfile,'a');
fprintf(fp,'filename\tt\tdiameter\tlink density\tself loops\tnum conn comp\tconnectivity\tsimple graph\tregular graph\tcomplete graph\teulerian graph\ttree\tgraphic\tmax node degree\tmin node degree\tmean node degree\tmax avg node degree\tmin avg node degree\tmean avg node degree\tclustering\tavg_degr\tmin(clsnss)\tmax(clsnss)\tmean(clsnss)\tmax(nb)\tmin(nb)\tmean(nb)\tmax(max((eb))\tmin(min((eb))\tmean(mean((eb))\tc2\tmax(wcc)\tmin(wcc)\tmean(wcc)\tapl\tsm_diam\tmax(vert_ecc)\tmin(vert_ecc)\tmean(vert_ecc)\tEntropyQ\tTheta\tSigma_0\tSigma_d\n');
%avg_degr,min(clsnss),max(clsnss),mean(clsnss),max(nb),max(max((eb)),min(min((eb)),mean(mean((eb)),c2,max(wcc),min(wcc),mean(wcc),apl,sm_diam,max(vert_ecc),min(vert_ecc),mean(vert_ecc)
fclose(fp);
for i=1:tstep:n-wind
    i
    n-wind
    y_fragm=y(i:i+wind);
    if(graph_type==1)
        Adj=double(crp(y_fragm,crp_emb_dimm,crp_delay,crp_epsilon));
    else
        Adj=ts2visgraph(y_fragm);
    end
    
    [d,dij]=diameter(Adj);%!!!!slow
    mas_diam=[mas_diam;d];
    ld=link_density(Adj); %easy
    mas_ld=[mas_ld;ld];
    sl=selfloops(Adj);
    mas_sl=[mas_sl;sl];
    ncc=num_conn_comp(Adj);%!
    mas_ncc=[mas_ncc;ncc];
    conn=isconnected(Adj);
    mas_conn=[mas_conn;conn];
    simplgr=issimple(Adj);
    mas_simplgr=[mas_simplgr;simplgr];
    reggr=isregular(Adj);
    mas_reggr=[mas_reggr;reggr];
    complgr=iscomplete(Adj);
    mas_complgr=[mas_complgr;complgr];
    eul_gr=iseulerian(Adj);
    mas_eulgr=[mas_eulgr;eul_gr];
    tree_gr=istree(Adj);
    mas_treegr=[mas_treegr;tree_gr];
    graphic_gr=isgraphic(Adj);
    mas_graphic_gr=[mas_graphic_gr;graphic_gr];
    %isbipartite(Adj)
    degr=degrees(Adj);
    mas_degr=[mas_degr;degr];
    maxdegr=max(degr);
    mas_maxdegr=[mas_maxdegr;maxdegr];
    mindegr=min(degr);
    mas_mindegr=[mas_mindegr;mindegr];
    meandegr=mean(degr);
    mas_meandegr=[mas_meandegr;meandegr];
    avrndegr=ave_neighbor_deg(Adj);
    mas_avrndegr=[mas_avrndegr;avrndegr];
    maxavrndegr=max(avrndegr);
    mas_maxavrndegr=[mas_maxavrndegr;maxavrndegr];
    minavrndegr=min(avrndegr);
    mas_minavrndegr=[mas_minavrndegr;minavrndegr];
    meanavrndegr=mean(avrndegr);
    mas_meanavrndegr=[mas_meanavrndegr;meanavrndegr];
    % edge_betweenness(Adj)% matrix
    clust=clust_coeff(Adj);
    mas_clust=[mas_clust;clust];
    % 20130418 Пропущенные меры
    % massivy:
    % mas_max_closeness,mas_min_node_betwenness,mas_max_edge_betwenness,mas_max_weighted_clust_coeff,mas_max_vertex_eccentricity
    avg_degr=average_degree(Adj);
    mas_average_degree=[mas_average_degree,avg_degr];
    clsnss=closeness(Adj,dij); %vector %!!!!!slow
    mas_max_closeness=[mas_max_closeness,max(clsnss)];
    mas_min_closeness=[mas_min_closeness,min(clsnss)];
    mas_mean_closeness=[mas_mean_closeness,mean(clsnss)];
    nb=node_betweenness_faster(Adj); %vector %!!!!!slow
    mas_max_node_betwenness=[mas_max_node_betwenness,max(nb)];
    mas_min_node_betwenness=[mas_min_node_betwenness,min(nb)];
    mas_mean_node_betwenness=[mas_mean_node_betwenness,mean(nb)];
    if (skip_complex_measures==0)
        eb=edge_betweenness(Adj); %matrix %!!!! slow
        mas_max_edge_betwenness=[mas_max_edge_betwenness,max(max(eb))];
        mas_min_edge_betwenness=[mas_min_edge_betwenness,min(min(eb))];
        mas_mean_edge_betwenness=[mas_mean_edge_betwenness,mean(mean(eb))];
    else
        %nb=zeros(length(Adj),1);
        %mas_max_node_betwenness=[mas_max_node_betwenness,0];
        %mas_min_node_betwenness=[mas_min_node_betwenness,0];
        %mas_mean_node_betwenness=[mas_mean_node_betwenness,0];
        eb=zeros(length(Adj),1);
        mas_max_edge_betwenness=[mas_max_edge_betwenness,0];
        mas_min_edge_betwenness=[mas_min_edge_betwenness,0];
        mas_mean_edge_betwenness=[mas_mean_edge_betwenness,0];
    end
    [c1,c2,c3]=clust_coeff(Adj);
    mas_clust2=[mas_clust2,c2];
    %mas_clust3=[mas_clust3,c3];
    wcc=weighted_clust_coeff(Adj); %vector
    mas_max_weighted_clust_coeff=[mas_max_weighted_clust_coeff,max(wcc)];
    mas_min_weighted_clust_coeff=[mas_min_weighted_clust_coeff,min(wcc)];
    mas_mean_weighted_clust_coeff=[mas_mean_weighted_clust_coeff,mean(wcc)];
    apl=ave_path_length(Adj,dij); %!!!slow
    mas_ave_path_length=[mas_ave_path_length,apl];
    sm_diam=smooth_diameter(Adj,thr_fr,dij); %!!!slow
    mas_smooth_diameter=[mas_smooth_diameter,sm_diam];
    vert_ecc=vertex_eccentricity(Adj,dij); % !!!!!slow ! max distance to any other node (vector)
    mas_max_vertex_eccentricity=[mas_max_vertex_eccentricity,max(vert_ecc)];
    mas_min_vertex_eccentricity=[mas_min_vertex_eccentricity,min(vert_ecc)];
    mas_mean_vertex_eccentricity=[mas_mean_vertex_eccentricity,mean(vert_ecc)];
    % avg_degr,min(clsnss),max(clsnss),mean(clsnss),max(nb),max(max((eb)),min(min((eb)),mean(mean((eb)),c2,max(wcc),min(wcc),mean(wcc),apl,sm_diam,max(vert_ecc),min(vert_ecc),mean(vert_ecc)
    % графика
    % draw_circ_graph
    % dot_matrix_plot()
    % radial_plot()
    % конец 20130418
    % начало 20130601 Энтропия сети
    if (skip_complex_measures==0)
        [theta,EntropyQ,W]=Theta(Adj,types,size(Adj,1),Q,Nrunmax);
        mas_theta=[mas_theta,theta];
        mas_EntropyQ=[mas_EntropyQ,EntropyQ];
        [Sigma_0,Sigma_d,W,z,kk]=Entropy_distance(Adj,dij,Nd);%,l);% ще є параметр l, який на даний момент задається в середині функції
        mas_Sigma_0=[mas_Sigma_0,Sigma_0];
        mas_Sigma_d=[mas_Sigma_d,Sigma_d];
    else
        theta=0;
        EntropyQ=0;
        Sigma_0=0;
        Sigma_d=0;
        mas_theta=[mas_theta,0];
        mas_EntropyQ=[mas_EntropyQ,0];
        mas_Sigma_0=[mas_Sigma_0,0];
        mas_Sigma_d=[mas_Sigma_d,0];
    end
    
    % Конец 20130601 Энтропия сети
    fp=fopen(logfile,'a');
    fprintf(fp,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',filename,i,d,ld,sl,ncc,conn,simplgr,reggr,complgr,eul_gr,tree_gr,graphic_gr,maxdegr,mindegr,meandegr,maxavrndegr,minavrndegr,meanavrndegr,clust,avg_degr,min(clsnss),max(clsnss),mean(clsnss),max(nb),min(nb),mean(nb),max(max(eb)),min(min(eb)),mean(mean(eb)),c2,max(wcc),min(wcc),mean(wcc),apl,sm_diam,max(vert_ecc),min(vert_ecc),mean(vert_ecc),EntropyQ,theta,Sigma_0,Sigma_d);
    fclose(fp);
    %link_density,selfloops,num_conn_comp ?find_conn_comp,giant_component,
    %?graph_complement,?tarjan,?isconnected,?issimple,
    %?isregular, ?iscomplete,?iseulerian, ?istree, ?isgraphic,
    %?isbipartite, add_edge_weights, degrees, ?laplacian_matrix,
    %ave_neighbor_deg, edge_betweenness, clust_coeff
    % 20130519 Добавляю локальные меры (сохранение каждого значения для
    % каждого узла)
    % degr - степень узла
    %clsnss - closeness
    %nb - node betweenness
    %eb - edge betweenness
    %c3??
    %wcc=weighted_clust_coeff
    %vert_ecc=vertex_eccentricity
    save_logs (logfile,i,degr,clsnss,nb,eb,c3,wcc,vert_ecc);
end

figure;
plot(mas_diam);
title('diameter');
outfile=strrep(filename,'.txt','_diameter.txt');
dlmwrite(outfile,mas_diam,'\n');


figure;
plot(mas_ld);
title('link density');
outfile=strrep(filename,'.txt','_link_density.txt');
dlmwrite(outfile,mas_ld,'\n');


%figure;
%plot(mas_sl);
%title('self loops number');
%outfile=strrep(filename,'.txt','_self_loops.txt');
%dlmwrite(outfile,mas_sl,'\n');

%figure;
%plot(mas_ncc);
%title('ncc');
%outfile=strrep(filename,'.txt','_num_conn_comp.txt');
%dlmwrite(outfile,mas_ncc,'\n');

%figure;
% plot(mas_conn);
% title('is the graph connected');
% outfile=strrep(filename,'.txt','_conn.txt');
% dlmwrite(outfile,mas_conn,'\n');

% figure;
% plot(mas_simplgr);
% title('is the graph simple');
% outfile=strrep(filename,'.txt','_simplegraph.txt');
% dlmwrite(outfile,mas_simplgr,'\n');
% 
% 
% 
% figure;
% plot(mas_complgr);
% title('is the graph complete');
% outfile=strrep(filename,'.txt','_completegraph.txt');
% dlmwrite(outfile,mas_complgr,'\n');
% 
% 
% figure;
% plot(mas_reggr);
% title('is the graph regular');
% outfile=strrep(filename,'.txt','_regular_graph.txt');
% dlmwrite(outfile,mas_reggr,'\n');
% 
% figure;
% plot(mas_eulgr);
% title('is the graph eulerian');
% outfile=strrep(filename,'.txt','_eulerian.txt');
% dlmwrite(outfile,mas_eulgr,'\n');
% 
% 
% figure;
% plot(mas_treegr);
% title('is the graph a tree');
% outfile=strrep(filename,'.txt','_tree.txt');
% dlmwrite(outfile,mas_treegr,'\n');
% 
% figure;
% plot(mas_graphic_gr);
% title('is the graph graphic');
% outfile=strrep(filename,'.txt','_graphic.txt');
% dlmwrite(outfile,mas_graphic_gr,'\n');

figure;
plot(mas_degr);
title('node degree dynamics');

figure;
plot(mas_maxdegr);
title('max node degree dynamics');
outfile=strrep(filename,'.txt','_maxdegr.txt');
dlmwrite(outfile,mas_maxdegr,'\n');


% figure;
% plot(mas_mindegr);
% title('min node degree dynamics');
% outfile=strrep(filename,'.txt','_mindegr.txt');
% dlmwrite(outfile,mas_mindegr,'\n');
% 
% 
% figure;
% plot(mas_meandegr);
% title('mean node degree dynamics');
% outfile=strrep(filename,'.txt','_meandegr.txt');
% dlmwrite(outfile,mas_meandegr,'\n');

figure;
plot(mas_avrndegr);
title('average node degree dynamics');

figure;
plot(mas_maxavrndegr);
title('max average node degree dynamics');
outfile=strrep(filename,'.txt','_maxavgndegr.txt');
dlmwrite(outfile,mas_maxavrndegr,'\n');


figure;
plot(mas_minavrndegr);
title('min average node degree dynamics');
outfile=strrep(filename,'.txt','_minavrndegr.txt');
dlmwrite(outfile,mas_minavrndegr,'\n');


figure;
plot(mas_meanavrndegr);
title('mean average node degree dynamics');
outfile=strrep(filename,'.txt','_meanavrndegr.txt');
dlmwrite(outfile,mas_meanavrndegr,'\n');


figure;
plot(mas_clust);
title('clustering coefficient');
outfile=strrep(filename,'.txt','_clustering.txt');
dlmwrite(outfile,mas_clust,'\n');

% figure;
% plot(mas_clust2);
% title('clustering coefficient2');
% outfile=strrep(filename,'.txt','_clustering2.txt');
% dlmwrite(outfile,mas_clust2,'\n');

figure;
plot(mas_average_degree);
title('Average degree');
outfile=strrep(filename,'.txt','_AvgDegree.txt');
dlmwrite(outfile,mas_average_degree,'\n');

% figure;
% plot(mas_max_closeness,'k');
% hold on;
% plot(mas_min_closeness,'k');
% plot(mas_mean_closeness,'r');
% title('closeness');
% outfile=strrep(filename,'.txt','_MaxCloseness.txt');
% dlmwrite(outfile,mas_max_closeness,'\n');
% outfile=strrep(filename,'.txt','_MinCloseness.txt');
% dlmwrite(outfile,mas_min_closeness,'\n');
% outfile=strrep(filename,'.txt','_MeanCloseness.txt');
% dlmwrite(outfile,mas_mean_closeness,'\n');

% figure;
% plot(mas_min_node_betwenness,'k');
% hold on;
% plot(mas_max_node_betwenness,'k');
% plot(mas_mean_node_betwenness,'r');
% title('Node betweenness');
% outfile=strrep(filename,'.txt','_MaxNodeBetweenness.txt');
% dlmwrite(outfile,mas_max_node_betwenness,'\n');
% outfile=strrep(filename,'.txt','_MinNodeBetweenness.txt');
% dlmwrite(outfile,mas_min_node_betwenness,'\n');
% outfile=strrep(filename,'.txt','_MeanNodeBetweenness.txt');
% dlmwrite(outfile,mas_mean_node_betwenness,'\n');

% figure;
% plot(mas_min_weighted_clust_coeff,'k');
% hold on;
% plot(mas_max_weighted_clust_coeff,'k');
% plot(mas_mean_weighted_clust_coeff,'r');
% title('Weighted Cluster coefficient');
% outfile=strrep(filename,'.txt','_MaxWeigthedClustCoef.txt');
% dlmwrite(outfile,mas_max_weighted_clust_coeff,'\n');
% outfile=strrep(filename,'.txt','_MinWeigthedClustCoef.txt');
% dlmwrite(outfile,mas_min_weighted_clust_coeff,'\n');
% outfile=strrep(filename,'.txt','_MeanWeigthedClustCoef.txt');
% dlmwrite(outfile,mas_mean_weighted_clust_coeff,'\n');

figure;
plot(mas_ave_path_length);
title('Average path length');
outfile=strrep(filename,'.txt','_AvgPathLength.txt');
dlmwrite(outfile,mas_ave_path_length,'\n');

% figure;
% plot(mas_smooth_diameter);
% title('Smooth Diameter');
% outfile=strrep(filename,'.txt','_SmoothDiameter.txt');
% dlmwrite(outfile,mas_smooth_diameter,'\n');

% figure;
% plot(mas_max_vertex_eccentricity,'k');
% hold on;
% plot(mas_min_vertex_eccentricity,'k');
% plot(mas_mean_vertex_eccentricity,'k');
% title('VertexExcentricity');
% outfile=strrep(filename,'.txt','_MaxVertexEccentricity.txt');
% dlmwrite(outfile,mas_max_vertex_eccentricity,'\n');
% outfile=strrep(filename,'.txt','_MinVertexEccentricity.txt');
% dlmwrite(outfile,mas_min_vertex_eccentricity,'\n');
% outfile=strrep(filename,'.txt','_MeanVertexEccentricity.txt');
% dlmwrite(outfile,mas_mean_vertex_eccentricity,'\n');

% figure;
% plot(mas_theta);
% title('\Theta');
% outfile=strrep(filename,'.txt','_Theta.txt');
% dlmwrite(outfile,mas_theta,'\n');
% 
% figure;
% plot(mas_EntropyQ);
% title('EntropyQ');
% outfile=strrep(filename,'.txt','_EntropyQ.txt');
% dlmwrite(outfile,mas_EntropyQ,'\n');
% 
% figure;
% plot(mas_Sigma_0);
% title('Sigma_0');
% outfile=strrep(filename,'.txt','_Sigma_0.txt');
% dlmwrite(outfile,mas_Sigma_0,'\n');
% 
% figure;
% plot(mas_Sigma_d);
% title('Sigma_d');
% outfile=strrep(filename,'.txt','_Sigma_d.txt');
% dlmwrite(outfile,mas_Sigma_d,'\n');

open_all_logs(logfile);
