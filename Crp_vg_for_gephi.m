% Graph from time series (CRP, Marwan)
infile='dax_91.txt';
outfile='dax_91_2000_3000_vizrgaph4gephi.csv';
type_graph=1; % 1 - crp, 2 - visibility graph
t_beg=2000;% Номер дня початку фрагмента
t_end=3000; % Номер дня кінця фрагмента
y=dlmread(infile);
n=length(y);
y_fragm=y(t_beg:t_end);
if (type_graph==1)
    X=double(crp(y_fragm));
else
    X=ts2visgraph(y_fragm);
end
[max_first,mas_last]=find(X==1);
linksdata=[max_first,mas_last];
ind_nonequal=find(max_first~=mas_last);
linksdata=linksdata(ind_nonequal,:);
[n1,m1]=size(linksdata);
linksdata2=[linksdata,ones(n1,1)];
dlmwrite(outfile,linksdata2,';');