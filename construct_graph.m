%construct_graph Constructs a graph of the dislocation network comprised of
%   edges which connect physical nodes
%   Arguments:
%   - links = links array containing node pairs which define each link
%   - conn = connectivity array constructed by the genconnectivity function
%   - physnodes = array of node numbers for the nodes which will terminate
%                   the edges
%
%   Outputs:
%   - G = consider deleting this!
%   - edges = similar to the links array but contains pairs of physical
%               nodes which bound graph edges
%   - edgenodes = array containing the set of all nodes contained in each
%               edge, each row is a separate edge, the first entry in each
%               row is the number of nodes in the edge
%   - edgelinks = same as edgenodes but gives the set of all links in
%               each edge

function [G,edges,edgenodes,edgelinks] = construct_graph(links,conn,physnodes)
linenum = 1;
edges = [];
nphys = length(physnodes);
for i=1:nphys
    nodei = physnodes(i);
    for j=1:conn(nodei,1)
        [linenodes,linelinks] = find_linenodes(nodei,j,links,conn);
        %only store the line if the starting physical node has a lower
        %index (to prevent duplicate storage)
        if linenodes(1)<linenodes(end)
            indend = find(linenodes(end)==physnodes);
            if isempty(indend)
                fprintf('Warning: isolated physical node detected\n')
                continue;
            end
            %edges(linenum,1:2) = [i indend];
            edges(linenum,1:2) = [linenodes(1) linenodes(end)];
            edgenodes(linenum,1:length(linenodes)+1) = [length(linenodes) linenodes];
            edgelinks(linenum,1:length(linelinks)+1) = [length(linelinks) linelinks];
            linenum = linenum+1;
        end
    end
end
%check for and remove duplicate edges
tmp = unique(edges,'rows');
if size(tmp,1)~=size(edges,1)
    fprintf('Warning: Duplicate links detected in construct_graph\n')
end
%edges = tmp;
%G = graph(edges(:,1),edges(:,2));
[~,G,~]=genconnectivity(zeros(max(links(:)),2),edges,10);
end

function [newnode,newlink] = find_nextnode(node,oldnode,branch,links,conne)
nlinks = conne(node,1);
onlinks = conne(node,2:2:2*nlinks);
if ismember(0,onlinks)
    onlinks
end
linknodes = links(onlinks,:);
%if an oldnode was provided, remove its links
if ~isempty(oldnode)
    %determine which link oldnode is on
    [oldlink,~] = find(linknodes==oldnode);
    %remove oldlink from the set of links
    linknodes(oldlink,:) = [];
    onlinks(oldlink) = [];
end
%grab nextnode and nextlink based on the specified branch
newnode = setdiff(linknodes(branch,:),node);
newlink = onlinks(branch);
end

function [linenodes,linelinks] = find_linenodes(node,arm,links,conn)
[newnode,newlink] = find_nextnode(node,[],arm,links,conn);
linenodes = [node newnode];
linelinks = [newlink];
while conn(newnode,1)==2
    oldnode = node;
    node = newnode;
    [newnode,newlink] = find_nextnode(node,oldnode,1,links,conn);
    linenodes = [linenodes newnode];
    linelinks = [linelinks newlink];
end
end