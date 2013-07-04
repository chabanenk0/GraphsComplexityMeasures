% Graph from time series (CRP, Marwan)
infile='dax_91.txt';
outfile='dax_91_2000_3000_vizrgaph4gephi.csv';
type_graph=1; % 1 - crp, 2 - visibility graph
t_beg=2000;% начало фрагмента
t_end=2050; % конец фрагмента
thr_fr=0.3;
crp_emb_dimm=1;
crp_delay=1;
crp_epsilon=0.5; %0.1 - Marwan default
y=dlmread(infile);
n=length(y);
y_fragm=y(t_beg:t_end);
if (type_graph==1) 
    X=double(crp(y_fragm,crp_emb_dimm,crp_delay,crp_epsilon));
else
    X=ts2visgraph(y_fragm);
end
i0=[];% вершина, которую выбирают как центр (пусто - по умолчанию)
figure;
if (isempty(i0))
    radial_plot(X)
else
    radial_plot(X,i0)
end
figure;
draw_circ_graph(X);

