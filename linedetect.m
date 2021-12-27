readnames=["init_ridges.png", "init_peaks.png", "init_valleys.png"; "post_ridges.png", "post_peaks.png", "post_valleys.png"];
savenames=["init_ridgenodes.mat", "init_peaknodes.mat", "init_valleynodes.mat"; "post_ridgenodes.mat", "post_peaknodes.mat", "post_valleynodes.mat"];
lseg=5;
lmin=5;

for ww=1:3
img = imread(readnames(qq,ww));
img = imcomplement(img);
if size(img,3)>1
    img = rgb2gray(img);
end
%figure(1)
%imshow(img)
%title('Original image')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Process and threshold the image
%convert to double for ease of processing
imgdbl = double(img);

%shift pixel values to be between 0 and 1
imin = min(min(imgdbl));
imax = max(max(imgdbl));
imgdbl = (imgdbl-imin)/(imax-imin);

    
    %flatten the image with the nlm filter
    imgdbl = imnlmfilt(imgdbl);

    %threshold the image using the adaptthresh function
    imgth = imgdbl;
    threshold = adaptthresh(imgth,0.5,'neighborhoodsize',2*floor(size(imgth)/64)+1);
    dinds = imgth>threshold;
    imgth(~dinds) = 0;
    imgth(dinds) = 1;

    %figure(2)
    %imshow(imgth)
    %title('Thresholded image')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract the skeleton
    %First we remove crap based on a minimum area criterion and that
    %circular blobs are removed
    imgskel = bwmorph(logical(imgth),'majority',Inf);
    imgskel = bwareafilt(imgskel,[100 Inf]);
    imgskel = bwpropfilt(imgskel,'extent',[0 0.5]);

    %first pass on constructing a skeleton
    imgskel = thin_pass(imgskel);


%figure(3);
tmp = imgdbl;
tmp(imgskel) = 1;
%h1 = imshowpair(imgskel,imgdbl,'diff');
%h1 = imshow(tmp);
%ax = gca;
%hold on

%Convert pixels to line segments
%one more thinning pass to smooth out issues with the manually drawn
%lines
%imgskel = thin_pass(imgskel);

%run a skeleton pass so we can use the endpoints and branchpoints
%commands reliably (as required by MATLAB)
imgskel = bwmorph(imgskel,'skel',Inf);

%identify branchpoint pixels
branchpoints = find(bwmorph(imgskel,'branchpoints'));
imgsize = size(imgskel);

%plot the final skeletonization
tmp = imgdbl;
tmp(imgskel) = 1;
%imshow(tmp);

%now loop over branchpoints and hop along branches to contstruct line
%segment representation of the network
[tmp1,tmp2] = ind2sub(imgsize,branchpoints);
rn = [tmp1 tmp2];
links = [];

%as we go, we remove pixels from the image so they don't get used again
imgcopy = imgskel;
for j=1:length(branchpoints)
    %grab the next branch point
    [I,J] = ind2sub(imgsize,branchpoints(j));

    %zero it out so we don't find it when be do neighbor searches
    imgcopy(I,J) = 0;
    hold on
    %plot(J,I,'rp','markersize',10)

    %loop until we explore all of the branchpoint's neighbors
    %we need to zero out the other neighbors to the branchpoint because
    %they can "pull" our search off of the intended branch using our
    %nearest neighbor algorithm
    
    %enumerate a list of all neighbors
    nbrlist = [];
    nbrcount = 0;
    while 1
        nbr = find_nbr(imgcopy,I,J);
        if isempty(nbr)
            break;
        end
        nbrcount = nbrcount+1;
        nbrlist(nbrcount,1:2) = nbr;
        %zero this neighbor out
        imgcopy(nbr(1),nbr(2)) = 0;
    end
    
    %go over all of the branches
    for k=1:nbrcount
        %turn the kth neighbor back on
        imgcopy(nbrlist(k,1),nbrlist(k,2)) = 1;
        [rn,links,imgcopy] = add_link(j,nbrlist(k,:),imgcopy,branchpoints,rn,links,lseg,lmin);
    end
    %while 1
    %    nbr = find_nbr(imgcopy,I,J);
    %    if isempty(nbr)
    %        break;
    %    end
    %    [rn,links,imgcopy] = add_link(j,nbr,imgcopy,branchpoints,rn,links,lseg,lmin);
    %end

    %toggle the branchpoint back on in case we encounter it in a
    %subsequent search
    imgcopy(I,J) = 1;
end
%all done with the branchpoints now, get rid of them
imgcopy(branchpoints) = 0;

%all we should have left now is isolated lines, which we can find based
%on their endpoints
endpoints = find(bwmorph(imgcopy,'endpoints'));
for j=1:length(endpoints)
    [I,J] = ind2sub(imgsize,endpoints(j));
    if imgcopy(I,J)==1
        imgcopy(I,J) = 0;
        rn = [rn; I J];
        nbr = find_nbr(imgcopy,I,J);
        if ~isempty(nbr)
            [rn,links,imgcopy] = add_link(size(rn,1),nbr,imgcopy,[],rn,links,lseg,lmin);
        end
    end
end

%now get rid of edges (connecting physical nodes) that are shorter than
%lmin
%loop over edges and determine their length
j = 1;
while 1
    %update the graph
    [conn,~,~]=genconnectivity(rn,links,10);
    physnodes = find(conn(:,1)~=2 & conn(:,1)>0); 
    [G,edges,edgenodes,edgelinks] = construct_graph(links,conn,physnodes);
    
    %check if we have tested all of the edges
    if j>size(edges,1)
        break;
    end
    
    %determine the length of this edge
    lenj = 0;
    Nnodesj = edgenodes(j,1);
    for k=1:Nnodesj-1
        lenj = lenj + norm(rn(edgenodes(j,k+1),:)-rn(edgenodes(j,k+2),:));
    end
    
    %if the edge is shorter than lmin, remove it by deleting all of its
    %links and merging the physical nodes together
    if lenj<lmin
        links(edgelinks(j,2:edgelinks(j,1)+1),:) = [];
        endnode1 = edgenodes(j,2);
        endnode2 = edgenodes(j,Nnodesj+1);
        %keep endnode1, goodbye to endnode2
        %move all of the connections to endnode1
        links(links==endnode2) = endnode1;
        %delete edgenode2, we have to renumber all of the nodes below it
        %by -1 
        belownodes = (links>endnode2);
        links(belownodes) = links(belownodes)-1;
        rn(endnode2,:) = [];
        
        %restart the search in case the order of the links was altered
        j = 1;
    else
        j = j+1;
    end
end

figure(4)
% imshow(imgdbl)
% hold on
% plotnodes2D([rn(:,2) rn(:,1)],links,'g');
% hold off

save(savenames(qq,ww),'rn','links')

end

function imgthin = thin_pass(img)
%first we thin it
imgthin = bwmorph(img,'thin',Inf);

%next fill and dilate to fill in any holes
imgthin = bwmorph(imgthin,'fill');
imgthin = imdilate(imgthin,[1 1 1; 1 1 1; 1 1 1]);

%thin again and clean it up (get rid of stray pixels)
imgthin = bwmorph(imgthin,'thin',Inf);
imgthin = bwmorph(imgthin,'clean');
end

%find a nearest neighbor pixel for pixel I J with preference toward 
%non-diagonal neighbors since these neighbors come first when traversing a
%skeleton
function nbr = find_nbr(img,I,J)
nbr = [];
%first check non-diagonals, then diagonals
is = [-1 1 0 0 -1 -1 1 1];
js = [0 0 -1 1 -1 1 -1 1];
for k=1:8
    i = is(k);
    j = js(k);
    if img(I+i,J+j)
        nbr = [I+i J+j];
        return;
    end
end
% for i=-1:1
%     for j=-1:1
%         if i==0 && j==0
%             continue;
%         end
%         if img(I+i,J+j)
%             nbr = [I+i J+j];
%             return;
%         end
%     end
% end
%now check diagonals
end

%add a link to the dislocation network starting at node1 and hopping along
%pixels starting with nbr
function [rn,links,img] = add_link(node1,nbr,img,branchpoints,rn,links,lseg,lmin)
imgsize = size(img);
count = 0;
len = 0;
Nn = size(rn,1);
rntmp = [];
linkstmp = [];

%zero out the current pixel
img(nbr(1),nbr(2)) = 0;
while 1
    oldnbr = nbr;
    nbr = find_nbr(img,nbr(1),nbr(2));
    
    %if we find an endpoint, add the endpoint node and the last segment
    if isempty(nbr)
        %only add a link if the total length exceeds lmin
        if len>=lmin
            rn = [rn; rntmp; oldnbr];
            links = [links; linkstmp; node1 Nn+1];
        end
        break;
    end
    
    %if we find another branchpoint, add the last segment using the
    %existing branchpoint node. we always add links which connect
    %branchpoints (no requirement on the minimum length)
    if ~isempty(branchpoints)
        tmp = ismember(branchpoints,sub2ind(imgsize,nbr(1),nbr(2)));
        if sum(tmp)==1
                rn = [rn; rntmp];
                links = [links; linkstmp; node1 find(tmp)];
            break;
        end
    end
    
    %this is just a pixel in the middle of the line, add a segment if
    %enough distance has passed since the last segment and then zero out
    %the pixel
    count = count+1;
    len = len+1;
    if count==lseg
        rntmp = [rntmp; nbr];
        Nn = Nn + 1;
        linkstmp = [linkstmp; node1 Nn];
        node1 = Nn;
        count = 0;
    end
    img(nbr(1),nbr(2)) = 0;
end
end

