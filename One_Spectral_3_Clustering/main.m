function [clusters3f,cheeger3f] = main(W)
% Computes a 3-partition of a given matrix by using eigenvectors of the 1-Laplacian.
%
% Input:
%   W: Sparse weight symmetric matrix.
%
% Output:
%   clusters3f: vector indicating the computed 3-clustering 
%   cheeger3f: the corresponding Cheeger constant
%
% (C)2019 Antonio Corbo Esposito, Domenico Angelo La Manna and Gianpaolo Piscitelli
% Dipartimento di Ingegneria Elettrica e dell'Informazione "M. Scarano",
% Via G. Di Biasio 43
% Università degli studi di Cassino e del Lazio Meridionale



n=size(W,1);
T=[0,0,0,0];
deg=full(sum(W,2));    
            normalized=false;
            [clusters2,~,cheegers2,v2] = OneSpectralClustering2nd(W,'ncc',2,0,0,2);
            [clusters3,~,cheegers3,v3] = OneSpectralClustering3rd(W,'ncc',2,0,0,2,v2);
            for i=1:n
                if (clusters2(i)==1)
                    if(clusters3(i)==1)
                        cluster(i)=1;
                        T(1)=1;
                    else
                        cluster(i)=2;
                        T(2)=2;
                    end
                else
                    if(clusters3(i)==1)
                        cluster(i)=3;
                        T(3)=3;
                    else
                        cluster(i)=4;
                        T(4)=4;
                    end
                end
            end          
            if(T>0)
                [cutpart1,cutpart2,cutpart3,cutpart4]=compute4thCutValue(cluster,W,normalized,deg);
                [~,index_max]=max([cutpart1,cutpart2,cutpart3,cutpart4]);
                cluster_new=[cluster' , cluster' , cluster'];
                for j=1:4
                    if (T(j)==index_max)
                        T(j)=0;
                        a=find(T>0);
                    end
                end
                for i=1:n
                    if (cluster(i)==index_max)
                        cluster_new(i,1)=a(1);
                        cluster_new(i,2)=a(2); 
                        cluster_new(i,3)=a(3);  
                    end
                end
                [cutpart11,cutpart12,cutpart13] = compute3rdCutValue(cluster_new(:,1),W,T,normalized,deg);                
                c1=max([cutpart11,cutpart12,cutpart13]);
                [cutpart21,cutpart22,cutpart23] = compute3rdCutValue(cluster_new(:,2),W,T,normalized,deg);
                c2=max([cutpart21,cutpart22,cutpart23]);
                [cutpart31,cutpart32,cutpart33] = compute3rdCutValue(cluster_new(:,3),W,T,normalized,deg);
                c3=max([cutpart31,cutpart32,cutpart33]);
                [cheeger3f,i_min]=min([c1,c2,c3]);
                clusters3f=cluster_new(:,i_min);
                i_miss=find(T==0);
                for i=1:n
                    if(clusters3f(i)==4)
                        clusters3f(i)=i_miss;
                    end
                end
            else
                R=find(T==0);
                S=size(R,2);
                if(S==1) 
                    clusters3f=cluster';
                    [cutpart1,cutpart2,cutpart3] = compute3rdCutValue(clusters3f,W,T,normalized,deg);
                    cheeger3f=max([cutpart1,cutpart2,cutpart3]);
                    i_miss=find(T==0);
                    for i=1:n
                        if(clusters3f(i)==4)
                            clusters3f(i)=i_miss;
                        end
                    end
                else
                   [clustersHB,~,cheegerHB] = OneSpectralClustering(W,'ncc',3,0,0,2); 
                   clusters3f=clustersHB(:,2);
%                   cheeger3f=cheegerHB(1,2);
T=[1 2 3];
[cutpartHB1,cutpartHB2,cutpartHB3] = compute3rdCutValue(clusters3f,W,T,false,deg);
cheeger3f=max([cutpartHB1,cutpartHB2,cutpartHB3]);
                end
            end   
end