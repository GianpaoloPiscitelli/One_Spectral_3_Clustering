function f=PseudoOrt(W,fold,f2,deg)
% Calculate a vector f (=u-lambda*f2, for suitable real lambda) such that 
% the pseudo scalar product between sign(f) and f2 is zero
%
% Input: 
%   fold: start vector
%   f2: the second eigenvector
%
% Output:
%   f: the vector that have zero pseudo product with f2
%
% (C)2020-21 Antonio Corbo Esposito and Gianpaolo Piscitelli
% Dipartimento di Ingegneria Elettrica e dell'Informazione "M. Scarano",
% Via G. Di Biasio 43
% Università degli studi di Cassino e del Lazio Meridionale
% https://github.com/GianpaoloPiscitelli/One_Spectral_3_Clustering


n=size(W);
t=0;
G=0;
I=0;
for i=1:n
    G=G+(abs(fold(i))*deg(i));
    for j=1:n
        I=I+W(i,j)*abs(fold(i)-fold(j));
    end
end
% t = starting value of the weighted-scalar product between sign(p) and v
for i = 1:n
    for j=1:n
       t=t+(W(i,j)*sign(fold(i)-fold(j))*(f2(i)-f2(j))*G);
    end
    t=t-(deg(i)*sign(fold(i))*f2(i)*I);%%%
end
%t=dot(sign(fold),f2)
if t<0
    f2=-f2;
    t=-t;
end

% computation of the third eigenvector candidate
if t==0 
    f=fold;
else
    %pr=fold.*f2 ;%%%
    for i=1:n    
        pr(i)=0;
    end
    for i = 1:n
        for j=1:n
            pr(i)=pr(i)+(W(i,j)*sign(fold(i)-fold(j))*(f2(i)-f2(j))*G);
        end
        pr(i)=pr(i)-(deg(i)*sign(fold(i))*f2(i)*I); %%%
    end
    % indices of positive products (concordances)
    index_c=find(pr>0); 
    v_c=f2(index_c);  
    pr_c=pr(index_c);
    % vector of concordances of u and v
    u_c=fold(index_c);  
    r=u_c./v_c;
    [r_ord,index]=sort(r);
    pr_p=pr_c(index);
    s=0;                      
    ii=1;
    % sum the ratios up to t/2
    while(s<t/2)
        s=s+pr_p(ii);
        %s=s+deg(ii)*abs(v_p(ii));%%%
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
    f=fold-lambda*f2;
end

%normalization of the third eigenvector candidate
%for i=1:n
%    G=G+(abs(f(i))*deg(i));
%end
%for i=1:n
%    f(i)=f(i)/G;
%end
end