function [Z, modes] = sr_zernike(x,y,K)
%SR_ZERNIKE Real, unit-disk RMS-normalized Zernike basis.
% Explicit ordering (NOT claimed to be Noll/Fringe): increasing radial n;
% within n, increasing |m|; m=0 first if present, then cos(m*t), sin(m*t).
% First four: 1, 2*x, 2*y, sqrt(3)*(2*(x.^2+y.^2)-1).
% Polynomial values outside disk are allowed; caller must mask apertures.
validateattributes(K,{'numeric'},{'scalar','integer','>=',1,'<=',153});
x=x(:); y=y(:); assert(numel(x)==numel(y),'Coordinate sizes differ.');
r=hypot(x,y); t=atan2(y,x); Z=zeros(numel(x),K); modes=zeros(K,3);
j=0; n=0;
while j<K
    for m=mod(n,2):2:n
        radial=zeros(size(r));
        for s=0:(n-m)/2
            c=(-1)^s*factorial(n-s)/(factorial(s)* ...
                factorial((n+m)/2-s)*factorial((n-m)/2-s));
            radial=radial+c*r.^(n-2*s);
        end
        if m==0
            j=j+1; Z(:,j)=sqrt(n+1)*radial; modes(j,:)=[n m 0];
        else
            j=j+1; Z(:,j)=sqrt(2*(n+1))*radial.*cos(m*t);
            modes(j,:)=[n m 1];
            if j<K
                j=j+1; Z(:,j)=sqrt(2*(n+1))*radial.*sin(m*t);
                modes(j,:)=[n m 2];
            end
        end
        if j>=K, break; end
    end
    n=n+1;
end
end
