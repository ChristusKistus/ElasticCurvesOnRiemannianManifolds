global G L K_0 T_0 Rtol;
Rtol=1e-5;
elastica=6;
mode=2

switch elastica
    case 0
        G=0; L=0; K_0=sqrt(20);
    case 1
        G=0; L=0; K_0=-10;
        % o=+-inf
    case 2
        G=-1; L=0; K_0=3;
        % -inf<o<-5.066
    case 2.1
        G=-1; L=0; K_0=sqrt(sqrt(40));
        % o=-6.32455 (-sqrt(40))
    case 3
        G=-2; L=0; K_0=-sqrt(-(G+L)*(5+1/15));
        % o=-5.066
    case 4
        G=-1*4; L=0; K_0=2.1*2;
        % -5.066<o<-4
    case 4.1
        G=-1; L=0; K_0=2.142;
        % o=-4.5882
    case 5
        G=-4; L=0; K_0=4;
        % o=-4
    case 6
        G=-1; L=0; K_0=1.7;
        % -4<o<-2
    case 7
        G=-2; L=0; K_0=2;
        % o=-2
    case 8
        G=-1; L=0; K_0=1;
        % -2<o<0
    case 9
        G=1; L=0; K_0=0;
        % o=0
    case 10
        G=1*25; L=0; K_0=1*5;
        % 0<o<inf
end
T_0=-diffcur(0)/sqrt(cur(0)^2+1)
% o=K_0^2/(G+L)
% z=[cx,cy,vx,vy,nx,ny] with c'=v %
% 0: test
% 1: right lintearia
% 2: squiggles
% 2.1: squiggles touching
% 3: cosine lemniscate
% 4: alternating squiggles
% 4.1: alternating squiggles touching
% 5: convict curve
% 6: pseudo trochoid facing down
% 7: circle
% 8: pseudo trochoid facing up
% 9: line
% 10: pseudo sinusoid

c_0=[0,0];
v_0=[1,1];
n_0=[-v_0(1,2),v_0(1,1)]*1/(sqrt(v_0(1,1).^2+v_0(1,2).^2));
n_dot0=-K_0*v_0;

C_0=[0,0,1];
V_0=[1,0,0];
N_0=[0,0,1];
B_0=[0,1,0];
N_dot0=-K_0*V_0+T_0*B_0;
B_dot0=-T_0*N_0;

switch mode
    case 1    
    initial1=[c_0(1,1),c_0(1,2),v_0(1,1),v_0(1,2),n_0(1,1),n_0(1,2)];
    initial2=[c_0(1,1),c_0(1,2),-v_0(1,1),-v_0(1,2),n_0(1,1),n_0(1,2)];
    opts = odeset('Reltol',1e-10,'AbsTol',1e-10,'Stats','on');
    tspan = [0 20];
    [t1,z1] = ode78(@fun, tspan, initial1, opts);
    [t2,z2] = ode78(@fun, tspan, initial2, opts);
    Curve=[flip(z1).',z2.'].';
    figure
    plot(Curve(:,1),Curve(:,2));
    title('elastica');
    case 2
    initial1=[1,0,0,1];
    initial2=[1,0,0,-1];
    %forman: [theta0,phi0,theta'0,phi'0]
    opts = odeset('Reltol',Rtol,'AbsTol',1e-10,'Stats','on');

    tspan = [0 10];
    [t1,y1] = ode78(@local1, tspan, initial1, opts);
    [t2,y2] = ode78(@local2, tspan, initial2, opts);
    T=[flip(-t1).',t2.'].';

    unprojcurve=([flip(y1).',y2.'].');
    Curve=zeros(size(unprojcurve,1),3);
    B=zeros(3,1);
    for i=1:size(unprojcurve,1)
        for j=1:3
            B=proj(unprojcurve(i,1),unprojcurve(i,2));
            Curve(i,j)=B(j,1);
        end
    end

    %figure
    %plot(T(:,1),sqrt(Curve(:,1).^2+Curve(:,2).^2+Curve(:,3).^2))

    figure
    sphere
    alpha 0.5, axis equal, view(3)
    hold on
    plot3(Curve(:,1),Curve(:,2),Curve(:,3),'LineWidth',4)
    title('3d elastica')

    case 3
    initial1=[1/sqrt(2),0,1/sqrt(2),0,1,0];
    initial2=[1/sqrt(2),0,1/sqrt(2),0,-1,0];
    opts = odeset('Reltol',1e-13,'AbsTol',1e-13,'Stats','on');

    tspan = [0 13];
    [t1,y1] = ode78(@guidobrunnet1, tspan, initial1, opts);
    [t2,y2] = ode78(@guidobrunnet2, tspan, initial2, opts);
    T=[flip(-t1).',t2.'].';

    Curve=[flip(y1).',y2.'].';
    %figure
    %plot(T(:,1),sqrt(Curve(:,1).^2+Curve(:,2).^2+Curve(:,3).^2))

    figure
    sphere
    alpha 0.5, axis equal, view(3)
    hold on
    plot3(Curve(:,1),Curve(:,2),Curve(:,3),'LineWidth',4);
    title('elastica');
    % distorted values for o. lemn: G=-0.76, K_0=sqrt(5.066)

end

function dz = fun(t,z)
dz = zeros(6,1);
dz(1) = z(3);
dz(2) = z(4);
dz(3) = cur(t)*z(5);
dz(4) = cur(t)*z(6);
dz(5) = -cur(t)*z(3);
dz(6) = -cur(t)*z(4);
end

% "Sphere" does not solve the problem. Somehow i got
% the torsion always wrong or smth. Doesn't matter
% anyway because the below version is better.
% w=[Cx,Cy,Cz,Vx,Vy,Vz,Nx,Ny,Nz,Bx,By,Bz]
% cur: geodesic curvature on whatever manifold with const gauß-curvature
% cur3d: actual curvature in R^3 on sphere
% tor: torsion of curve on sphere

function dw = Sphere(t,w)
dw = zeros(12,1); 
dw(1) = w(4);
dw(2) = w(5);
dw(3) = w(6);
dw(4) = cur3d(t)*w(7);
dw(5) = cur3d(t)*w(8);
dw(6) = cur3d(t)*w(9);
dw(7) = -cur3d(t)*w(4)+tor(t)*w(10);
dw(8) = -cur3d(t)*w(5)+tor(t)*w(11);
dw(9) = -cur3d(t)*w(6)+tor(t)*w(12);
dw(10) = -tor(t)*w(7);
dw(11) = -tor(t)*w(8);
dw(12) = -tor(t)*w(9);
end

% some random diff eq for curves on spheres i found
% in a research paper by guido brunnet 
% (works pretty much just as fine as local)

function dy=guidobrunnet1(t,y)
dy=zeros(6,1);
dy(1)=y(4);
dy(2)=y(5);
dy(3)=y(6);
dy(4)=-y(1)-cur(t)*(y(2)*y(6)-y(3)*y(5));
dy(5)=-y(2)-cur(t)*(y(3)*y(4)-y(1)*y(6));
dy(6)=-y(3)-cur(t)*(y(1)*y(5)-y(2)*y(4));
end

function dy=guidobrunnet2(t,y)
dy=zeros(6,1);
dy(1)=y(4);
dy(2)=y(5);
dy(3)=y(6);
dy(4)=-y(1)+cur(t)*(y(2)*y(6)-y(3)*y(5));
dy(5)=-y(2)+cur(t)*(y(3)*y(4)-y(1)*y(6));
dy(6)=-y(3)+cur(t)*(y(1)*y(5)-y(2)*y(4));
end

%relations: theta''*f1+phi''*g1+r1=0
%           theta''*f2+phi''*g2+r2=0,
% "local" because this solves the curve in local
% coordinates: theta,phi

function dy=local1(t,y)
f1=@(x,y) -sin(x)*cos(y);
g1=@(x,y) -cos(x)*sin(y);
f2=@(x,y) -sin(x)*sin(y);
g2=@(x,y) cos(x)*cos(y);
r1=@(x,y,xp,yp,t) sin(x)*cos(x)*yp^2*f1(x,y)-2*tan(x)*xp*yp*g1(x,y)+cur(t)*(xp*sin(y)-yp*sin(x)*cos(x)*cos(y))/sqrt(xp^2+yp^2*cos(x)^2);
r2=@(x,y,xp,yp,t) sin(x)*cos(x)*yp^2*f2(x,y)-2*tan(x)*xp*yp*g2(x,y)+cur(t)*(-xp*cos(y)-yp*sin(x)*cos(x)*sin(y))/sqrt(xp^2+yp^2*cos(x)^2);

dy = zeros(4,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),t)-g2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(4) = (f1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),t)-f2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2)));
end

function dy=local2(t,y)
f1=@(x,y) -sin(x)*cos(y);
g1=@(x,y) -cos(x)*sin(y);
f2=@(x,y) -sin(x)*sin(y);
g2=@(x,y) cos(x)*cos(y);
r1=@(x,y,xp,yp,t) sin(x)*cos(x)*yp^2*f1(x,y)-2*tan(x)*xp*yp*g1(x,y)-cur(t)*(xp*sin(y)-yp*sin(x)*cos(x)*cos(y))/sqrt(xp^2+yp^2*cos(x)^2);
r2=@(x,y,xp,yp,t) sin(x)*cos(x)*yp^2*f2(x,y)-2*tan(x)*xp*yp*g2(x,y)-cur(t)*(-xp*cos(y)-yp*sin(x)*cos(x)*sin(y))/sqrt(xp^2+yp^2*cos(x)^2);

dy = zeros(4,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),t)-g2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(4) = (f1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),t)-f2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2)));
end

% y=[theta1,phi1,theta2,phi2], yp=[theta1',phi1',theta2',phi2'] aka [theta1',phi1',theta1'',phi1''] with theta1,phi1 being the actualy coordinates, and
% theta1'=theta2, phi1'=phi2

% res equations: (yp(3)+sin(y(1))*cos(y(1)))*(-sin(y(1))*cos(y(2)))+(yp(4)-2*tan(y(1))*yp(1)*yp(2))*(-cos(y(1))*sin(y(2)))-cur(t)*(cos(y(1))*cos(y(2)))
%        (yp(3)+sin(y(1))*cos(y(1)))*(-sin(y(1))*sin(y(2)))+(yp(4)-2*tan(y(1))*yp(1)*yp(2))*(cos(y(1))*sin(y(2)))-cur(t)*(cos(y(1))*sin(y(2)))
%        (yp(3)+sin(y(1))*cos(y(1)))*(cos(y(1)))+(yp(4)-2*tan(y(1))*yp(1)*yp(2))*(sin(y(1)))-cur(t)*sin(y(1))];

function PROJ=proj(theta,phi)
    PROJ=[cos(theta)*cos(phi),cos(theta)*sin(phi),sin(theta)].';
end

function CUR3D=cur3d(s)
    CUR3D=sqrt(cur(s)^2+1);
end

function DIFFCUR=diffcur(s)
    global Rtol;
    DIFFCUR=(cur(s+Rtol)-cur(s))/Rtol;
end

function TOR=tor(s)
    TOR=diffcur(s)/(1+sqrt(cur(s)^2+1));
end

function CUR=cur(s)
global G L K_0;

if G+L>=0
    r=sqrt(G+L+(K_0.^2)/2);
    k=abs(K_0)/(2*r);
    CUR=sign(K_0)*K_0*cn(r*s,k);
elseif G+L<0
    if K_0==abs(sqrt(-2*(G+L)))
        CUR=sqrt(-2*(G+L));
    elseif K_0==abs(2*sqrt(-(G+L)))
        CUR=(2*sign(K_0)*sqrt(-(G+L)))/cosh(s*sqrt(-(G+L)));
    elseif abs(K_0)<sqrt(-2*(G+L))
        r=sqrt(-4*(G+L)-K_0.^2)/2;
        k=sqrt(2*r^2+G+L)/r;
        f=@(x) 1./sqrt(1-k.^2*sin(x).^2);
        s0=integral(f,0,pi/2)/r;
        CUR=sign(K_0)*2*r*dn(r*(s+s0),k);
    elseif sqrt(-2*(G+L))<abs(K_0) && abs(K_0)<2*sqrt(-(G+L))
        r=abs(K_0)/2;
        k=sqrt(2*r^2+G+L)/r;
        CUR=sign(K_0)*2*r*dn(r*s,k);
    elseif abs(K_0)>2*sqrt(-(G+L))
        r=sqrt(G+L+(K_0.^2)/2);
        k=abs(K_0)/(2*r);
        CUR=sign(K_0)*K_0*cn(r*s,k);
    end
end
end

function CN=cn(s,k)
[~,C,~]=ellipj(s,k^2);
CN=C;
end

function DN=dn(s,k)
[~,~,D]=ellipj(s,k^2);
DN=D;

end
