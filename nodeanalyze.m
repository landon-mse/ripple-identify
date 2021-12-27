loadnames=["init_ridgenodes.mat", "init_peaknodes.mat", "init_valleynodes.mat"; "post_ridgenodes.mat", "post_peaknodes.mat", "post_valleynodes.mat"];
for n=1:3
load(loadnames(qq,n))
xlength=0;
ylength=0; 
nodes = 1:size(rn,1);
LINKMAX=length(links(:,1));
for q=1:LINKMAX
    n0=links(q,1);
    n1=links(q,2);
    if ismember(n0,nodes) && ismember(n1,nodes)
        ydist = abs(rn(n0,2)-rn(n1,2));
        xdist = abs(rn(n0,1)-rn(n1,1));
        xlength = xlength + xdist;
        ylength = ylength + ydist;
    end
end

myvalues(n,1) = xlength; 
myvalues(n,2) = ylength;

end