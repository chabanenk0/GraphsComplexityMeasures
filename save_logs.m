function save_logs (baselogfile,i,degr,clsnss,nb,eb,c3,wcc,vert_ecc);
% Функция сохранения (дополнения) логов для локальных мер (тех, что выдают
% характеристику для каждого значения.
% i - 
% degr - степень узла
% clsnss - closeness
% nb - node betweenness
% eb - edge betweenness
% c3??
% wcc=weighted_clust_coeff
% vert_ecc=vertex_eccentricity

n=length(degr);
logfile=strrep(baselogfile,'.txt','_degrees.txt');
save_log_series(logfile,i,degr);
logfile=strrep(baselogfile,'.txt','_nb.txt');
save_log_series(logfile,i,nb);
logfile=strrep(baselogfile,'.txt','_eb.txt');
save_log_series(logfile,i,eb);
logfile=strrep(baselogfile,'.txt','_c3.txt');
save_log_series(logfile,i,c3);
logfile=strrep(baselogfile,'.txt','_wcc.txt');
save_log_series(logfile,i,wcc);
logfile=strrep(baselogfile,'.txt','_vert_ecc.txt');
save_log_series(logfile,i,vert_ecc);
