function [y,z]=phi_radau(R,l,mu,t)
D=[sqrt(mu)*R;eye(l+1)];
d=[zeros(l,1);eye(l+1,1)];
%z=inv(mu*R'*R+eye(l+1,l+1))*eye(l+1,1);
[Q,R]=qr(D);
z=R\(Q'*d);
y=t*z'*z;
end