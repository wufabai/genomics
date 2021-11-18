%% This one goes through All the clusters that has ERPs in them and count their origins
clear all;
cd('');
load('');
load('Asg8Cluster24x80homcut.mat');
A=Asg8Cluster;
B=asg8selfblast;
Raln1=B.length./B.slen;
Raln2=B.length./B.qlen;
Raln3=max(Raln1,Raln2);
Rscore=Raln3.*B.pident;
B.Rscore=Rscore;

MAGlist=[]; % presence in each MAG
MAGbase=[0 0 0 0 0 0 0 0];
PHYlist=[]; % presence in each phylum
PHYbase=[0 0 0 0 0];
EUKCLcomp=[]; % In the clusters that contain ERPs who else are there
EUKCLbase=[0 0 0 0];
% find clusters that contain Euk
ind1=find(A.TAXSN==1);
Eukcl1=unique(A.CLSN(ind1));
TabEukExtended=table();
% analyse cluster by cluster
for i=1:numel(Eukcl1)
    ind2=find(A.CLSN==Eukcl1(i));
    Tabcl1=A(ind2,:);
    % find nonEuk ones
    ind3=find(Tabcl1.TAXSN>1);
    ind4=find(Tabcl1.TAXSN==0);
    
    if isempty(ind3) && ~isempty(ind4) % only with unknown
        Tabcl1.TAXSN(ind4)=1;
    end
    
    if ~isempty(ind3) && ~isempty(ind4)
        % Now we have to check to which one they are the closest
        for j=1:numel(ind4)
            header1=string(Tabcl1.Header(ind4(j)));
            ind5=find(strcmp(B.qseqid,header1));
            scores=[];
            for k=1:numel(ind5)
                header2=string(B.sseqid(ind5(k)));
                Rscore2=B.Rscore(ind5(k));
                ind6=find(strcmp(A.Header,header2));
                TAXSN2=A.TAXSN(ind6);
                scores=cat(1,scores,[Rscore2 TAXSN2]);
            end
            sortscores=sortrows(scores,1,'descend');
            ind7=find(sortscores(:,2)>0);
            if isempty(ind7)
                Tabcl1.TAXSN(ind4(j))=0;
            else
                Tabcl1.TAXSN(ind4(j))=sortscores(ind7(1),2);
            end
        end
    end
        %% now assign to MAGs/Phylum
        EUKCLbase1=EUKCLbase;
        PHYbase1=PHYbase;
        MAGbase1=MAGbase;
    for m=1:numel(ind2)
        taxon=Tabcl1.TAXSN(m);
        EUKCLbase1(taxon+1)=EUKCLbase1(taxon+1)+1;
        if taxon==1
           MAGbase1(Tabcl1.MAGSN(m))=MAGbase1(Tabcl1.MAGSN(m))+1;
           PHYbase1(Tabcl1.PHYSN(m))=PHYbase1(Tabcl1.PHYSN(m))+1;
        end
    end
    MAGlist=cat(1,MAGlist,MAGbase1);
    PHYlist=cat(1,PHYlist,PHYbase1);
    EUKCLcomp=cat(1,EUKCLcomp,EUKCLbase1);
    
    indeuk2=find(Tabcl1.TAXSN==1);
    TabEukExtended=cat(1,TabEukExtended,Tabcl1(indeuk2,:));
end
PHYlistb=PHYlist;
PHYlistb(PHYlistb>1)=1;    
PHYsum=sum(PHYlistb,2);
histogram(PHYsum);
ind8=find(PHYsum==5);
ind9=(Eukcl1(ind8));
AA=TabEukExtended;
AAA=sortrows(AA,21);
writetable(AAA,'EUKclusters_20210329_24x80homcut.xlsx');    
    
    
    
    