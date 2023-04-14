function [y,s]=dphi_gauss(R,l,z,mu,t)

%s=inv(mu*R'*R+eye(l,l))*z;
D=[sqrt(mu)*R;eye(l)];
d=[zeros(l,1);z];
[Q,R]=qr(D);
s=R\(Q'*d);
y=(2*t/mu)*z'*(s-z);
end