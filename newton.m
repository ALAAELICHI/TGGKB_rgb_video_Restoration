function mu=newton(R,l,t,mu,eta,s)
 while phi_gauss(R,l,mu,t)>s^2+1/10*(eta^2-1)*s^2 
   [y1,z]=phi_gauss(R,l,mu,t);
    [y2,w]=dphi_gauss(R,l,z,mu,t);
    y3=ddphi_gauss(R,l,z,w,mu,t);
    betta=mu+(y2/(2*y3));
alpha=2*y2*sqrt(mu-betta);
gamma=y1-alpha*sqrt(mu-betta);
mu=betta+gamma^2/alpha^2;
end     
end