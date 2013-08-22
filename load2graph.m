function res = load2graph(filename)

y = dlmread(filename);

n = length(y);
wind = 500;
tstep = 5;
thr_fr = 0.5;
crp_emb_dimm = 1;
crp_delay = 1;
crp_epsilon = 0.4; %0.1 - Marwan default
graph_type = 2; %1 - crp, 2 - visibility

Q = 5;
types = round((1:wind+1)/(wind/(Q-1))+1);
Nrunmax = 3;
Nd = 10;
skip_complex_measures = 1;

mas_diam = [];
mas_degr = [];
mas_maxdegr = [];
mas_mindegr = [];
mas_meandegr = [];
mas_Hurst = [];
mas_MFDFA_width = [];
mas_MFDFA_max = [];

logfilename=strrep(filename,'.txt','_logfile.txt');
fp=fopen(logfilename,'a');
fprintf(fp,'t\tmas_diam\tmas_maxdegr\tmas_mindegr\tmas_meandegr\tmas_Hurst\tmas_MFDFA_width\tmas_MFDFA_max\n');
fclose(fp);


h = waitbar(0,'Calculating graphs');
for i = 1:tstep:n-wind % главный цикл
    waitbar(i/(n-wind),h);
    y_fragm = y(i:i+wind);
    if(graph_type==1)
        Adj=double(crp(y_fragm,crp_emb_dimm,crp_delay,crp_epsilon));
    else
        Adj=ts2visgraph(y_fragm);
    end
    
    [d,dij]=diameter(Adj);
    mas_diam=[mas_diam;d];

    degr = degrees(Adj);
    mas_degr = [mas_degr;degr];
    maxdegr = max(degr);
    mas_maxdegr = [mas_maxdegr; maxdegr];
    mindegr = min(degr);
    mas_mindegr = [mas_mindegr;mindegr];
    meandegr = mean(degr);
    mas_meandegr = [mas_meandegr;meandegr];
    % ----- Hurst
    [hx,hy] = Hurst2(degr,50,1.1);
    hx = log(hx); hy = log(hy);
    [a,b] = MinSquare(hx,hy);
    mas_Hurst = [mas_Hurst; a];
    % ----- MFDFA
%    calculateMFDFA(degr,20,400,[-3:0.1:-0.1 0.1:0.1:3],2,'nRatio','100','90');
    [S,FqS,Q]=MFDFA(degr,20,400,[-3:0.1:-0.1 0.1:0.1:3],2,'nRatio','100');
    [Hq,TAUq,Q2]=MFDFAHTau(S,FqS,Q);
    [XResTemp,YResTemp]=MFDFAAlphaFTAUSpline(Q2,TAUq,'90');
    mas_MFDFA_width = [mas_MFDFA_width; (max(XResTemp)-min(XResTemp))];
    [maxMFDFA,imaxMFDFA] = max(YResTemp);
    mas_MFDFA_max = [mas_MFDFA_max; XResTemp(imaxMFDFA(1))];
    clsnss=closeness(Adj,dij); %vector %!!!!!slow
    nb=node_betweenness_faster(Adj); %vector %!!!!!slow
    if (skip_complex_measures==0)
        eb=edge_betweenness(Adj); %matrix %!!!! slow
    else
        %nb=zeros(length(Adj),1);
        %mas_max_node_betwenness=[mas_max_node_betwenness,0];
        %mas_min_node_betwenness=[mas_min_node_betwenness,0];
        %mas_mean_node_betwenness=[mas_mean_node_betwenness,0];
        eb=zeros(length(Adj),1);
    end
    [c1,c2,c3]=clust_coeff(Adj);
    wcc=weighted_clust_coeff(Adj); %vector
    vert_ecc=vertex_eccentricity(Adj,dij); % !!!!!slow ! max distance to any other node (vector)
    fp=fopen(logfilename,'a');
    fprintf(fp,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',i,mas_diam,mas_maxdegr,mas_mindegr,mas_meandegr,mas_Hurst,mas_MFDFA_width,mas_MFDFA_max);
    fclose(fp);
    save_logs (logfilename,i,degr,clsnss,nb,eb,c3,wcc,vert_ecc);
end
close(h);

res = cell(7,2);
res{1,1} = mas_diam;      % диаметр графа
res{1,2} = 'DIAMETER';
res{2,1} = mas_maxdegr;   % максимальные степени вершин
res{2,2} = 'MAXDEGR ';
res{3,1} = mas_mindegr;   % минимальные степени вершин
res{3,2} = 'MINDEGR ';
res{4,1} = mas_meandegr;  % среднее значение степеней вершин
res{4,2} = 'MEANDEGR';
res{5,1} = mas_Hurst;
res{5,2} = 'HURST   ';
res{6,1} = mas_MFDFA_width;
res{6,2} = 'DFAWIDTH';
res{7,1} = mas_MFDFA_max;
res{7,2} = 'DFAMAX  ';

figure;
plot(mas_diam);
title('diameter');
outfile=strrep(filename,'.txt','_diameter.txt');
dlmwrite(outfile,mas_diam,'\n');

figure;
plot(mas_degr);
title('node degree dynamics');

figure;
plot(mas_maxdegr);
title('max node degree dynamics');
outfile=strrep(filename,'.txt','_maxdegr.txt');
dlmwrite(outfile,mas_maxdegr,'\n');

figure;
plot(mas_maxdegr);
title('degree Hurst');
outfile=strrep(filename,'.txt','_maxdegr.txt');
dlmwrite(outfile,mas_maxdegr,'\n');

figure;
plot(mas_Hurst);
title('degree Hurst');
outfile=strrep(filename,'.txt','_DegreeHurst.txt');
dlmwrite(outfile,mas_Hurst,'\n');

figure;
plot(mas_MFDFA_width);
title('degree Hurst');
outfile=strrep(filename,'.txt','_DegreeMFDFA_width.txt');
dlmwrite(outfile,mas_MFDFA_width,'\n');

figure;
plot(mas_MFDFA_max);
title('degree Hurst');
outfile=strrep(filename,'.txt','_DegreeMFDFA_max.txt');
dlmwrite(outfile,mas_MFDFA_max,'\n');

open_all_logs(logfile);
