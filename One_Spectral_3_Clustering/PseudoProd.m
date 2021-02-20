function f=PseudoProd(W,fold,eigv,deg)
% Calculate a vector f (=u-lambda*v, for suitable real lambda) such that 
% the weighted scalar product between sign(f) and v is zero
%
% Input: 
%   fold: start vector
%   eigv: the second eigenvector
%
% Output:
%   f: the vector that have the weighted-null scalar product of his sign with v
%
% (C)2020-21 Antonio Corbo Esposito and Gianpaolo Piscitelli
% Dipartimento di Ingegneria Elettrica e dell'Informazione "M. Scarano",
% Via G. Di Biasio 43
% Università degli studi di Cassino e del Lazio Meridionale
% https://github.com/GianpaoloPiscitelli/One_Spectral_3_Clustering



n=size(W);
d=sum(W);
t=0;
G=0;
% t = starting value of the weighted-scalar product between sign(p) and v
%for i = 1:n
%    t=t+(d(i)*sign(fold(i))*eigv(i));
%end
t=dot(sign(fold),eigv);
if t<0
    eigv=-eigv;
    t=-t;
end

% computation of the third eigenvector candidate
if t==0 
    f=fold;
else
    pr=fold.*eigv ;
    % indices of positive products (concordances of u and v)
    index_c=find(pr>0); 
    v_c=eigv(index_c);   
    % vector of concordances of u and v
    u_c=fold(index_c);   
    r=u_c./v_c;
    [r_ord,index]=sort(r);
    v_p=v_c(index);
    s=0;                      
    ii=1;
    % sum the ratios up to t/2
    while(s<t/2)
        s=s+abs(v_p(ii));
        ii=ii+1;
    end
    if (s>t/2)
        pos=ii-1;
        pos2=ii-1;
    else
        pos=ii-1;
        pos2=ii;
    end
    % estimate the parameter lambda
    lambda=(r_ord(pos)+r_ord(pos2))/2;
    % calculate the hypotized function
    f=fold-lambda*eigv;
end

%normalization of the third eigenvector candidate
%for i=1:n
%    G=G+(abs(f(i))*deg(i));
%end
%for i=1:n
%    f(i)=f(i)/G;
%end
end