%% A greedy clustering approach for all proteins in the 8 selected Asgard archaea genomes.
%% Read fasta containing all proteins and record serial number and their MAG assignment (1-8)
%% Read the master table from Diamond blast, and do an identity and length filter
%% Save all data into 1 table
%% list follow faa list, column 1: MAG SN. 2. Cluster SN
% clear all;
cd('');
load('');

homcut=24; % homology score threshold for clustering
Lcut=0.8; % length threshold
Scut=100; % must be 100a.a. long

A=asgard8annotation;
% alignment length filtered
B=asg8selfblast;
Raln1=B.length./B.slen;
Raln2=B.length./B.qlen;
Raln3=max(Raln1,Raln2);
Rscore=Raln3.*B.pident;
B.Rscore=Rscore;
Raln4=intersect(find(Rscore>=homcut),find(Raln3>=Lcut));
Raln5=intersect(find(B.slen>Scut),find(B.qlen>Scut));

Raln=intersect(Raln4,Raln5);
Tab=B(Raln,:);

% load and log the headers
Nheaders=size(A,1);
A.PEPSN=(1:1:Nheaders)';
A.MAGSN=zeros(Nheaders,1);
A.PHYSN=zeros(Nheaders,1);
for i=1:size(MAGassign,1)
    MAGstr=string(MAGassign.MAGname(i));
    ind0=find(contains(asgard8annotation.Header,MAGstr));
    A.MAGSN(ind0)=MAGassign.MAGSN(i);
    A.PHYSN(ind0)=MAGassign.PHYLUMSN(i);
end

% Start making clusters
cl=1;
MAGSNtemp=A.MAGSN;
MAGSNlog=[MAGSNtemp MAGSNtemp];
MAGSNlog(:,2)=0;
pos=1;
Na=0;
Header2=A.Header(pos);
while sum(MAGSNtemp)>0
    
    HeaderCluster=Header2;
    ind1=find(strcmp(Tab.qseqid,Header2));
    if ~isempty(ind1)
        HeaderCluster=cat(1,HeaderCluster, string(Tab.sseqid(ind1)));
    end
    ind2=find(strcmp(Tab.sseqid,Header2));
    if ~isempty(ind2)
        HeaderCluster=cat(1,HeaderCluster, string(Tab.qseqid(ind2)));
    end
    HeaderClusterb=unique(HeaderCluster);
    Nb=numel(HeaderClusterb);
    %% Now go through the the list find new ones associated with the cluster
    %% Until no new ones get recruited.
    while Nb>Na        
        for j=Na+1:Nb
            Header2=HeaderClusterb(j);
            ind1=find(strcmp(Tab.qseqid,Header2));
            if ~isempty(ind1)
                HeaderCluster=cat(1,HeaderCluster, string(Tab.sseqid(ind1)));
            end
            ind2=find(strcmp(Tab.sseqid,Header2));
            if ~isempty(ind2)
                HeaderCluster=cat(1,HeaderCluster, string(Tab.qseqid(ind2)));
            end
        end
        HeaderClusterb=unique(HeaderCluster);
        Na=Nb;
        Nb=numel(HeaderClusterb);
    end
    %% Now go recruit the serial numbers and deplete them from the list.
    for k=1:Nb
        ind3=find(strcmp(A.Header,HeaderClusterb(k)));
        MAGSNtemp(ind3)=0;
        MAGSNlog(ind3,2)=cl;
    end
    cl = cl+1; % next cluster
    ind4=find(MAGSNtemp);
    if ~isempty(ind4)
        Header2=A.Header(ind4(1));
    end
end
A.CLSN=MAGSNlog(:,2);

%% Finally, I serialize all the taxonomic assignments
%% 0, undefined, 1, euk, 2, bac, 3, arc
A.TAXSN=ones(Nheaders,1);
indbac=find(ismember(A.TaxonScope,'Bacteria'));
indarc=find(ismember(A.TaxonScope,'Archaea'));
indunk=find(isundefined(A.TaxonScope));
A.TAXSN(indbac)=2;
A.TAXSN(indarc)=3;
A.TAXSN(indunk)=0;
Asg8Cluster=A;    
save('Asg8Cluster24x80homcut.mat','Asg8Cluster');
            
