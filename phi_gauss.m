function [y,z]=phi_gauss(R,l,mu,t)
 %z=inv(mu*R'*R+eye(l,l))*eye(l,1);
 D=[sqrt(mu)*R;eye(l)];
d=[zeros(l,1);eye(l,1)];
[Q,R]=qr(D);
z=R\(Q'*d);
y=t*z'*z;
end