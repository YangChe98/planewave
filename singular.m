clc
clear
 k1=1;
 k2=2;
 k3=3;
  
%  intf8= integral3(f,-pi,0,-pi,0,-pi,0);
%   intf5= integral3(f,-pi,0,0,pi,0,pi);
%   intf6= integral3(f,-pi,0,0,pi,-pi,0);
%   intf7= integral3(f,-pi,0,-pi,0,0,pi);
%  g=@(x1,y1,z1) f(pi*x1,pi*y1,pi*z1)*pi^3;
%   intg=integral3(g,0,1,0,1,0,1);
%   h=@(x,u,w) g(x,x.*u,x.*w);
%   syms x2 y2 z2
%   g(x2,y2,z2)
intf=[];
%error=[];
K=40;
% for k1=-K/2:K/2-1
%     for k2=-K/2:K/2-1
%         for k3=-K/2:K/2-1
            for k1=0:K
    for k2=0:K
        for k3=0:K
            k1 
            k2
            k3
       %     f=@(x,y,z) exp(1i*(k1*x+k2*y+k3*z))./(x.^2+y.^2+z.^2).^(1/2);
  %intf1=integral3(f,-pi,pi,-pi,pi,-pi,pi)
%   f1=@(x,y,z) exp(1i*(k1*x*pi+k2*y*pi+k3*z*pi));
%   g1=@(x,u,w) (f1(x,x.*u,x.*w)+ f1(x.*u,x.*w,x)+ f1(x.*w,x,x.*u))*pi^2;
% h1=@(x,u,w) g1(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth1=integral3(h1,0,1,0,1,0,1);
% 
%  f2=@(x,y,z) exp(1i*(k1*x*pi+k2*y*pi-k3*z*pi));
%   g2=@(x,u,w) (f2(x,x.*u,x.*w)+ f2(x.*u,x.*w,x)+ f2(x.*w,x,x.*u))*pi^2;
% h2=@(x,u,w) g2(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth2=integral3(h2,0,1,0,1,0,1);
% 
%  f3=@(x,y,z) exp(1i*(k1*x*pi-k2*y*pi+k3*z*pi));
%   g3=@(x,u,w) (f3(x,x.*u,x.*w)+ f3(x.*u,x.*w,x)+ f3(x.*w,x,x.*u))*pi^2;
% h3=@(x,u,w) g3(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth3=integral3(h3,0,1,0,1,0,1);
% 
%  f4=@(x,y,z) exp(1i*(k1*x*pi-k2*y*pi-k3*z*pi));
%  g4=@(x,u,w) (f4(x,x.*u,x.*w)+ f4(x.*u,x.*w,x)+ f4(x.*w,x,x.*u))*pi^2;
% h4=@(x,u,w) g4(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth4=integral3(h4,0,1,0,1,0,1);
% 
%   f5=@(x,y,z) exp(1i*(-k1*x*pi+k2*y*pi+k3*z*pi));
%   g5=@(x,u,w) (f5(x,x.*u,x.*w)+ f5(x.*u,x.*w,x)+ f5(x.*w,x,x.*u))*pi^2;
% h5=@(x,u,w) g5(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth5=integral3(h5,0,1,0,1,0,1);
% 
%  f6=@(x,y,z) exp(1i*(-k1*x*pi+k2*y*pi-k3*z*pi));
%   g6=@(x,u,w) (f6(x,x.*u,x.*w)+ f6(x.*u,x.*w,x)+ f6(x.*w,x,x.*u))*pi^2;
% h6=@(x,u,w) g6(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth6=integral3(h6,0,1,0,1,0,1);
% 
%  f7=@(x,y,z) exp(1i*(-k1*x*pi-k2*y*pi+k3*z*pi));
%   g7=@(x,u,w) (f7(x,x.*u,x.*w)+ f7(x.*u,x.*w,x)+ f7(x.*w,x,x.*u))*pi^2;
% h7=@(x,u,w) g7(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth7=integral3(h7,0,1,0,1,0,1);
% 
%  f8=@(x,y,z) exp(1i*(-k1*x*pi-k2*y*pi-k3*z*pi));
%  g8=@(x,u,w) (f8(x,x.*u,x.*w)+ f8(x.*u,x.*w,x)+ f8(x.*w,x,x.*u))*pi^2;
% h8=@(x,u,w) g8(x,u,w).*x./(1+u.^2+w.^2).^(1/2);
%  inth8=integral3(h8,0,1,0,1,0,1);
%  inth=inth1+inth2+inth3+inth4+inth5+inth6+inth7+inth8
h1=@(x,u,w) pi^2*(exp(1i*(k1*x*pi+k2*x.*u*pi+k3*x.*w*pi))+exp(1i*(k1*x.*u*pi+k2*x.*w*pi+k3*x*pi))+exp(1i*(k1*x.*w*pi+k2*x*pi+k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth1=integral3(h1,0,1,0,1,0,1);


h2=@(x,u,w) pi^2*(exp(1i*(k1*x*pi+k2*x.*u*pi-k3*x.*w*pi))+exp(1i*(k1*x.*u*pi+k2*x.*w*pi-k3*x*pi))+exp(1i*(k1*x.*w*pi+k2*x*pi-k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth2=integral3(h2,0,1,0,1,0,1);


h3=@(x,u,w) pi^2*(exp(1i*(k1*x*pi-k2*x.*u*pi+k3*x.*w*pi))+exp(1i*(k1*x.*u*pi-k2*x.*w*pi+k3*x*pi))+exp(1i*(k1*x.*w*pi-k2*x*pi+k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth3=integral3(h3,0,1,0,1,0,1);


h4=@(x,u,w) pi^2*(exp(1i*(k1*x*pi-k2*x.*u*pi-k3*x.*w*pi))+exp(1i*(k1*x.*u*pi-k2*x.*w*pi-k3*x*pi))+exp(1i*(k1*x.*w*pi-k2*x*pi-k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth4=integral3(h4,0,1,0,1,0,1);

h5=@(x,u,w) pi^2*(exp(1i*(-k1*x*pi+k2*x.*u*pi+k3*x.*w*pi))+exp(1i*(-k1*x.*u*pi+k2*x.*w*pi+k3*x*pi))+exp(1i*(-k1*x.*w*pi+k2*x*pi+k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth5=integral3(h5,0,1,0,1,0,1);

 
h6=@(x,u,w) pi^2*(exp(1i*(-k1*x*pi+k2*x.*u*pi-k3*x.*w*pi))+exp(1i*(-k1*x.*u*pi+k2*x.*w*pi-k3*x*pi))+exp(1i*(-k1*x.*w*pi+k2*x*pi-k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth6=integral3(h6,0,1,0,1,0,1);

 
h7=@(x,u,w) pi^2*(exp(1i*(-k1*x*pi-k2*x.*u*pi+k3*x.*w*pi))+exp(1i*(-k1*x.*u*pi-k2*x.*w*pi+k3*x*pi))+exp(1i*(-k1*x.*w*pi-k2*x*pi+k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth7=integral3(h7,0,1,0,1,0,1);

h8=@(x,u,w) pi^2*(exp(1i*(-k1*x*pi-k2*x.*u*pi-k3*x.*w*pi))+exp(1i*(-k1*x.*u*pi-k2*x.*w*pi-k3*x*pi))+exp(1i*(-k1*x.*w*pi-k2*x*pi-k3*x.*u*pi))).*x./(1+u.^2+w.^2).^(1/2);
 inth8=integral3(h8,0,1,0,1,0,1);
 inth=inth1+inth2+inth3+inth4+inth5+inth6+inth7+inth8
       intf=[intf;inth];
      % error=[error;intf1-inth];
        end
    end
end
