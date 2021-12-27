function p0=plotnodes2D(rn,links,color,nodes)
%plot nodes

if(nargin < 4)
    nodes = 1:size(rn,1);
end

plot(0,0); hold on;
LINKMAX=length(links(:,1));
for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    if ismember(n0,nodes) && ismember(n1,nodes)
        p0=plot(rn([n0,n1],1),rn([n0,n1],2),strcat(color,'.-'));
    end
end
hold off
axis equal