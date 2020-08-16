function f=FalseProd(W,u,v)
% Calculate a vector f (=u-lambda*v, for suitable real lambda) such that 
% the weighted scalar product between sign(f) and v is zero
%
% Input: 
%   u: start vector
%   v: the vector such that f will have the null weighted-scalar product (between
%      sign(f) and v)
%
% Output:
%   f: the vector that have the weighted-null scalar product of his sign with v
%
% (C)2019 Antonio Corbo Esposito, Domenico Angelo La Manna and Gianpaolo Piscitelli
% Dipartimento di Ingegneria Elettrica e dell'Informazione "M. Scarano",
% Via G. Di Biasio 43
% Università degli studi di Cassino e del Lazio Meridionale


% t = starting value of the weighted-scalar product between sign(p) and v
n=size(W);
d=sum(W);
t=0;
for (i = 1:n)
    t=t+(d(i)*u(i)*v(i));
end
t=dot(sign(u),v);
if t<0
    v=-v;
    t=-t;
end


if t==0 
    f=u;
else
    pr=u.*v ;
    % indices of positive products (concordances of u and v)
    index_c=find(pr>0); 
    v_c=v(index_c);   
    % vector of concordances of u and v
    u_c=u(index_c);   
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
    f=u-lambda*v; 
    %n=sum(abs(f))
end
end