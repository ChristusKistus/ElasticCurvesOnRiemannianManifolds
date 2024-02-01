global G L K_0 T_0 Rtol Dtol P sur;
Rtol=1e-5;
Dtol=1e-5;
elastica=5;

% \ONLY PLANAR!\
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

mode=5;
% 1: planar, 
% 2: sphere (using local chart),
% 3: sphere (guidobrunnets ode in global coords),
% 4: general surface with G=0,
% 5: general surface.
sur=10;

% 1: sphere
% 2: cylinder
% 3: hyperboloid (1 shell)
% 4: möbius
% 5: "klein bottle"
% 6: cos rotation
% 7: torus
% 8: hyperbolic paraboloid
% 9: pseudosphere (tractroid)
% 10: screw
% 11: chocolate
cameraflight=0;

switch elastica
    case 0
        G=0; L=0; K_0=sqrt(20);
    case 1
        G=0; L=0; K_0=-10;
        % o=+-inf
    case 2
        G=0; L=-1; K_0=3;
        % -inf<o<-5.066
    case 2.1
        G=0; L=-1; K_0=sqrt(sqrt(40));
        % o=-6.32455 (-sqrt(40))
    case 3
        G=0; L=-2; K_0=-sqrt(-(G+L)*(5+1/15));
        % o=-5.066
    case 4
        G=0; L=-1*4; K_0=2.1*2;
        % -5.066<o<-4
    case 4.1
        G=0; L=-1; K_0=2.142;
        % o=-4.5882
    case 5
        G=0; L=-4; K_0=4;
        % o=-4
    case 6
        G=0; L=-1; K_0=1.7;
        % -4<o<-2
    case 7
        G=0; L=-2; K_0=2;
        % o=-2
    case 8
        G=0; L=-1; K_0=1;
        % -2<o<0
    case 9
        G=0; L=1; K_0=0;
        % o=0
    case 10
        G=0; L=1*25; K_0=1*5;
        % 0<o<inf
end
T_0=-diffcur(0)/sqrt(cur(0)^2+1)

c_0=[0,0];
v_0=[1,0];
n_0=[-v_0(1,2),v_0(1,1)]*1/(sqrt(v_0(1,1).^2+v_0(1,2).^2));
n_dot0=-K_0*v_0;

C_0=[0,0,1];
V_0=[1,0,0];
N_0=[0,0,1];
B_0=[0,1,0];
N_dot0=-K_0*V_0+T_0*B_0;
B_dot0=-T_0*N_0;

switch mode
    case 4
    case 5 
        switch sur
            case 0
                initial1=[0.1,0,0,1,1,0,K_0,0];
                initial2=[0.1,0,0,1,-1,0,K_0,0];
            case 1
                initial1=[0.1,0,0,1,1,0,K_0,0];
                initial2=[0.1,0,0,-1,-1,0,K_0,0];
            case 2
                initial1=[0,0,1,0,0,1,K_0,0];
                initial2=[0,0,-1,0,0,-1,-K_0,0];
            case 3
                initial1=[0.1,0,0,1,1,0,K_0,0];
                initial2=[0.1,0,0,-1,-1,0,-K_0,0];
            case 4
                initial1=[0.1,0,0,1,1,0,K_0,0];
                initial2=[0.1,0,0,-1,-1,0,-K_0,0];
            case 5
                initial1=[0.1,0,0,1,1,0,K_0,0];
                initial2=[0.1,0,0,-1,-1,0,-K_0,0];
            case 6
                initial1=[0.3,0,0,1,1,0,K_0,0];
                initial2=[0.3,0,0,-1,-1,0,-K_0,0];
            case 7
                initial1=[0,pi/2,0,1,1,0,K_0,0];
                initial2=[0,pi/2,0,-1,-1,0,-K_0,0];
            case 8
                initial1=[0,0,0,1,1,0,K_0,0];
                initial2=[0,0,0,-1,-1,0,-K_0,0];
            case 9
                initial1=[0,2,0,1,1,0,K_0,0];
                initial2=[0,2,0,-1,-1,0,-K_0,0];            
            case 10
                initial1=[0.5,1,0,1,1,0,K_0,0];
                initial2=[0.5,1,0,-1,-1,0,-K_0,0];            
            case 11
                initial1=[0.5,1,0,1,1,0,K_0,0];
                initial2=[0.5,1,0,-1,-1,0,-K_0,0];
        end
end
%}

switch mode
    case 1    
    initial11=[c_0(1,1),c_0(1,2),v_0(1,1),v_0(1,2),n_0(1,1),n_0(1,2)];
    initial12=[c_0(1,1),c_0(1,2),-v_0(1,1),-v_0(1,2),n_0(1,1),n_0(1,2)];
    opts = odeset('Reltol',1e-10,'AbsTol',1e-10,'Stats','on');
    tspan = [0 20];
    [t1,z1] = ode78(@fun, tspan, initial11, opts);
    [t2,z2] = ode78(@fun, tspan, initial12, opts);
    Curve=[flip(z1).',z2.'].';
    figure
    plot(Curve(:,1),Curve(:,2));
    title('elastica');
    case 2
    initial1=[0.1,0,0,1,1,0];
    initial2=[0.1,0,0,-1,-1,0];
    %format: [theta0,phi0,theta'0,phi'0,eta1,eta2]
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
    sphere, alpha 0.5, axis equal, view(3)
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
    case 4
    initial1=[0,0,1,0,0,1];
    initial2=[0,0,-1,0,0,-1];
    opts = odeset('Reltol',Rtol,'AbsTol',1e-5,'Stats','on');

    v=diffU(initial1(1),initial1(2));
    [MAX,I]=max([abs([1 0 0]*cross(v(:,1),v(:,2))) 
                 abs([0 1 0]*cross(v(:,1),v(:,2))) 
                 abs([0 0 1]*cross(v(:,1),v(:,2)))]);
    switch I
        case 1
            P=[2;3];
        case 2
            P=[1;3];
        case 3
            P=[1;2];
    end
    tspan = [0 5];
    [t1,y1] = ode78(@manifold1, tspan, initial1, opts);
    [t2,y2] = ode78(@manifold2, tspan, initial2, opts);
    T=[flip(-t1).',t2.'].';

    unprojcurve=([flip(y1).',y2.'].');
    Curve=zeros(size(unprojcurve,1),3);
    B=zeros(3,1);
    for i=1:size(unprojcurve,1)
        for j=1:3
            B=param(unprojcurve(i,1),unprojcurve(i,2));
            Curve(i,j)=B(j,1);
        end
    end

    %figure
    %plot(T(:,1),sqrt(Curve(:,1).^2+Curve(:,2).^2+Curve(:,3).^2))

    figure
    plot3(Curve(:,1),Curve(:,2),Curve(:,3),'LineWidth',4),axis equal, view(3)
    title('3d elastica')
    hold on
    M = 30; N = 30;
    xpar = linspace(-pi,pi,M); 
    ypar = linspace(-1,1,N); 
    [XPAR YPAR] = meshgrid(xpar,ypar);
    X=zeros(size(XPAR,1),size(YPAR,1));    
    Y=zeros(size(XPAR,1),size(YPAR,1));
    Z=zeros(size(XPAR,1),size(YPAR,1));
    for i=1:size(XPAR)
        for j=1:size(YPAR)
            Coord=param(XPAR(i,j),YPAR(i,j));
            X(i,j)=Coord(1);
            Y(i,j)=Coord(2);
            Z(i,j)=Coord(3);
        end
    end
    surf(X,Y,Z), alpha 0.3

    case 5
    % initial: [x0,y0,x0',y0',eta1,eta2,K0,K0']
    % x,y being the local coordinates, eta1, eta2 the spanning vectors for
    % the normal vector, K0,K0' curvature initial values
    % !!!!!as such, the initial values x0,y0 always have to be perpendicular
    % to eta1,eta2 with eta1=-eta2!!!!!
    

    v=diffU(initial1(1),initial1(2));
    [MAX,I]=max([abs([1 0 0]*cross(v(:,1),v(:,2))) 
             abs([0 1 0]*cross(v(:,1),v(:,2))) 
             abs([0 0 1]*cross(v(:,1),v(:,2)))]);
    switch I
        case 1
            P=[2;3];
        case 2
            P=[1;3];
        case 3
            P=[1;2];
    end

    opts = odeset('Reltol',Rtol,'AbsTol',1e-5,'Stats','on');

    tspan = [0 14];
    [t1,y1] = ode78(@surface, tspan, initial1, opts);
    [t2,y2] = ode78(@surface, tspan, initial2, opts);
    T=[flip(-t1).',t2.'].';

    unprojcurve=([flip(y1).',y2.'].');
    Curve=zeros(size(unprojcurve,1),3);
    B=zeros(3,1);
    for i=1:size(unprojcurve,1)
        for j=1:3
            B=param(unprojcurve(i,1),unprojcurve(i,2));
            Curve(i,j)=B(j,1);
        end
    end

    %figure
    %plot(T(:,1),sqrt(Curve(:,1).^2+Curve(:,2).^2+Curve(:,3).^2))

    figure
    plot3(Curve(:,1),Curve(:,2),Curve(:,3),'LineWidth',4),axis equal, view(3)
    title('3d elastica')
    hold on
    M = 50; N = 50;
    xpar = linspace(-pi,pi,M); 
    ypar = linspace(-pi,pi,N); 
    [XPAR YPAR] = meshgrid(xpar,ypar);
    X=zeros(size(XPAR,1),size(YPAR,1));    
    Y=zeros(size(XPAR,1),size(YPAR,1));
    Z=zeros(size(XPAR,1),size(YPAR,1));
    for i=1:size(XPAR)
        for j=1:size(YPAR)
            Coord=param(XPAR(i,j),YPAR(i,j));
            X(i,j)=Coord(1);
            Y(i,j)=Coord(2);
            Z(i,j)=Coord(3);
        end
    end
    D0=zeros(M,N,3);
    for i=1:M
        for j=1:N
            U=param(xpar(i),ypar(j));
            D0(i,j,1)=30/360;
            D0(i,j,2)=1;
            D0(i,j,3)=0.2+(4+U(3))/10;
        end

    end
    C0=hsv2rgb(D0);
    switch sur
        case 11
            S=surf(X,Y,Z,C0);
            S.FaceAlpha=0.8;
        otherwise
            surf(X,Y,Z), alpha 0.3
    end
end
% camera flight
switch cameraflight
    case 1
    switch mode
        case {2,3,4,5}
        camproj perspective
        camva(55)
        theta_dat=1/3*sin(linspace(0,pi,100))+0.2;
        phi_dat=linspace(0,2*pi,100);
        X_dat=4*sur*cos(theta_dat).*sin(phi_dat);
        Y_dat=4*sur*cos(theta_dat).*cos(phi_dat);
        Z_dat=4*sur*sin(theta_dat);
        U_dat=-X_dat;
        V_dat=-Y_dat;
        W_dat=-Z_dat;
        for i=1:length(X_dat)
            campos([X_dat(i),Y_dat(i),Z_dat(i)]);
            camtarget([U_dat(i),V_dat(i),W_dat(i)]);
            drawnow;
        end
    end
end

function U=param(x,y)
global sur
U=zeros(3,1);

switch sur
    case 1
        U(1)=cos(x)*cos(y);
        U(2)=cos(x)*sin(y);
        U(3)=sin(x);
    case 2
        U(1)=cos(x);
        U(2)=sin(x);
        U(3)=y;
    case 3
        U(1)=3*cos(x)*cosh(y);
        U(2)=3*sin(x)*cosh(y);
        U(3)=3*sinh(y);
    case 4
        U(1)=0.8*cos(x)*(1+y/2*cos(x/2));
        U(2)=0.8*sin(x)*(1+y/2*cos(x/2));
        U(3)=0.8*y/2*sin(x/2);
    case 5
        if x<=pi
            U(1)=0.5*(3*(1+sin(x)+2*(1-cos(x)/2)*cos(y)))*cos(x);
            U(2)=0.5*(4+2*(1-cos(x)/2)*cos(y))*sin(x);
            U(3)=0.5*2*(1-cos(x)/2)*sin(y);
        else
            U(1)=0.5*3*(1+sin(x))*cos(x)-2*(1-cos(x)/2)*cos(y);
            U(2)=0.5*4*sin(x);
            U(3)=0.5*2*(1-cos(x)/2)*sin(y);
        end
    case 6
        U(1)=0.8*x;
        U(2)=0.8*y;
        U(3)=0.4*cos(3*sqrt(x^2+y^2));
    case 7
        U(1)=(1+1/2*cos(y))*cos(x);
        U(2)=(1+1/2*cos(y))*sin(x);
        U(3)=1/2*sin(y);
    case 8
        U(1)=x;
        U(2)=y;
        U(3)=(x^2-y^2);
    case 9
        U(1)=20*(sech(y))*cos(x);
        U(2)=20*(sech(y))*sin(x);
        U(3)=20*(y-tanh(y));
    case 10
        U(1)=3*x*cos(y);
        U(2)=3*x*sin(y);
        U(3)=3*y;
    case 11
        U(1)=4*x;
        U(2)=4*y;
        U(3)=4*(10*sin(2*x)/(abs(10*sin(2*x))+1))*(10*sin(2*y)/(abs(10*sin(2*y))+1));
end
%}
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
% y=[theta,phi,theta',phi',eta1,eta2] with eta1, eta2 spanning the normal
% vector: n=eta1*dF/dtheta+eta2*dF/dphi

function dy=local1(t,y)
f1=@(x,y) -sin(x)*cos(y);
g1=@(x,y) -cos(x)*sin(y);
f2=@(x,y) -sin(x)*sin(y);
g2=@(x,y) cos(x)*cos(y);
r1=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*yp^2*f1(x,y)-2*tan(x)*xp*yp*g1(x,y)+cur(t)*(-eta1*sin(x)*cos(y)-eta2*cos(x)*sin(y));
r2=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*yp^2*f2(x,y)-2*tan(x)*xp*yp*g2(x,y)+cur(t)*(-eta1*sin(x)*sin(y)+eta2*cos(x)*cos(y));
r3=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*eta2*yp*f1(x,y)-tan(x)*(eta1*yp+eta2*xp)*g1(x,y)-cur(t)*(-xp*sin(x)*cos(y)-yp*cos(x)*sin(y));
r4=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*eta2*yp*f2(x,y)-tan(x)*(eta1*yp+eta2*xp)*g2(x,y)-cur(t)*(-xp*sin(x)*sin(y)+yp*cos(x)*cos(y));

dy = zeros(6,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),y(5),y(6),t)-g2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),y(5),y(6),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(4) = (f1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),y(5),y(6),t)-f2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),y(5),y(6),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2)));
dy(5) = (g1(y(1),y(2))*r4(y(1),y(2),y(3),y(4),y(5),y(6),t)-g2(y(1),y(2))*r3(y(1),y(2),y(3),y(4),y(5),y(6),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(6) = (f1(y(1),y(2))*r4(y(1),y(2),y(3),y(4),y(5),y(6),t)-f2(y(1),y(2))*r3(y(1),y(2),y(3),y(4),y(5),y(6),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2))); 
end

function dy=local2(t,y)
f1=@(x,y) -sin(x)*cos(y);
g1=@(x,y) -cos(x)*sin(y);
f2=@(x,y) -sin(x)*sin(y);
g2=@(x,y) cos(x)*cos(y);
r1=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*yp^2*f1(x,y)-2*tan(x)*xp*yp*g1(x,y)-cur(t)*(-eta1*sin(x)*cos(y)-eta2*cos(x)*sin(y));
r2=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*yp^2*f2(x,y)-2*tan(x)*xp*yp*g2(x,y)-cur(t)*(-eta1*sin(x)*sin(y)+eta2*cos(x)*cos(y));
r3=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*eta2*yp*f1(x,y)-tan(x)*(eta1*yp+eta2*xp)*g1(x,y)+cur(t)*(-xp*sin(x)*cos(y)-yp*cos(x)*sin(y));
r4=@(x,y,xp,yp,eta1,eta2,t) sin(x)*cos(x)*eta2*yp*f2(x,y)-tan(x)*(eta1*yp+eta2*xp)*g2(x,y)+cur(t)*(-xp*sin(x)*sin(y)+yp*cos(x)*cos(y));

dy = zeros(6,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),y(5),y(6),t)-g2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),y(5),y(6),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(4) = (f1(y(1),y(2))*r2(y(1),y(2),y(3),y(4),y(5),y(6),t)-f2(y(1),y(2))*r1(y(1),y(2),y(3),y(4),y(5),y(6),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2)));
dy(5) = (g1(y(1),y(2))*r4(y(1),y(2),y(3),y(4),y(5),y(6),t)-g2(y(1),y(2))*r3(y(1),y(2),y(3),y(4),y(5),y(6),t))/(g2(y(1),y(2))*f1(y(1),y(2))-g1(y(1),y(2))*f2(y(1),y(2)));
dy(6) = (f1(y(1),y(2))*r4(y(1),y(2),y(3),y(4),y(5),y(6),t)-f2(y(1),y(2))*r3(y(1),y(2),y(3),y(4),y(5),y(6),t))/(f2(y(1),y(2))*g1(y(1),y(2))-f1(y(1),y(2))*g2(y(1),y(2))); 
end

function dy=manifold1(t,y)
global P
DF=diffU(y(1),y(2));
CS=christoffel(y(1),y(2));

f1=DF(P(1),1);
f2=DF(P(2),1);
g1=DF(P(1),2);
g2=DF(P(2),2);
r1=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f1+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g1-cur(t)*(y(5)*f1+y(6)*g1);
r2=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f2+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g2-cur(t)*(y(5)*f2+y(6)*g2);
r3=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f1+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g1+cur(t)*(y(3)*f1+y(4)*g1);
r4=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f2+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g2+cur(t)*(y(3)*f2+y(4)*g2);

dy = zeros(6,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1*r2-g2*r1)/(g2*f1-g1*f2);
dy(4) = (f1*r2-f2*r1)/(f2*g1-f1*g2);
dy(5) = (g1*r4-g2*r3)/(g2*f1-g1*f2);
dy(6) = (f1*r4-f2*r3)/(f2*g1-f1*g2);
end


function dy=manifold2(t,y)
global P
DF=diffU(y(1),y(2));
CS=christoffel(y(1),y(2));

f1=DF(P(1),1);
f2=DF(P(2),1);
g1=DF(P(1),2);
g2=DF(P(2),2);
r1=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f1+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g1+cur(t)*(y(5)*f1+y(6)*g1);
r2=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f2+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g2+cur(t)*(y(5)*f2+y(6)*g2);
r3=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f1+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g1-cur(t)*(y(3)*f1+y(4)*g1);
r4=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f2+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g2-cur(t)*(y(3)*f2+y(4)*g2);

dy = zeros(6,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1*r2-g2*r1)/(g2*f1-g1*f2);
dy(4) = (f1*r2-f2*r1)/(f2*g1-f1*g2);
dy(5) = (g1*r4-g2*r3)/(g2*f1-g1*f2);
dy(6) = (f1*r4-f2*r3)/(f2*g1-f1*g2);
end

% y=[u1,u2,u1',u2',eta1,eta2,K1,K2]
function dy=surface(t,y)
global L P;
DF=diffU(y(1),y(2));
CS=christoffel(y(1),y(2));

f1=DF(P(1),1);
f2=DF(P(2),1);
g1=DF(P(1),2);
g2=DF(P(2),2);
r1=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f1+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g1-y(7)*(y(5)*f1+y(6)*g1);
r2=(CS(1,1)*y(3)*y(3)+CS(1,2)*y(3)*y(4)+CS(2,1)*y(4)*y(3)+CS(2,2)*y(4)*y(4))*f2+(CS(1,3)*y(3)*y(3)+CS(1,4)*y(3)*y(4)+CS(2,3)*y(4)*y(3)+CS(2,4)*y(4)*y(4))*g2-y(7)*(y(5)*f2+y(6)*g2);
r3=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f1+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g1+y(7)*(y(3)*f1+y(4)*g1);
r4=(CS(1,1)*y(5)*y(3)+CS(1,2)*y(5)*y(4)+CS(2,1)*y(6)*y(3)+CS(2,2)*y(6)*y(4))*f2+(CS(1,3)*y(5)*y(3)+CS(1,4)*y(5)*y(4)+CS(2,3)*y(6)*y(3)+CS(2,4)*y(6)*y(4))*g2+y(7)*(y(3)*f2+y(4)*g2);

dy = zeros(8,1);
dy(1) = y(3);
dy(2) = y(4);
dy(3) = (g1*r2-g2*r1)/(g2*f1-g1*f2);
dy(4) = (f1*r2-f2*r1)/(f2*g1-f1*g2);
dy(5) = (g1*r4-g2*r3)/(g2*f1-g1*f2);
dy(6) = (f1*r4-f2*r3)/(f2*g1-f1*g2);
dy(7) = y(8);
dy(8) = -1/2*y(7)^3-(GC(y(1),y(2))+L)*y(7);
end

function D=diffU(x,y)
Dtol=1e-8;
    D1=(param(x+Dtol,y)-param(x,y))/Dtol;
    D2=(param(x,y+Dtol)-param(x,y))/Dtol;
    D=[D1,D2];
end

% diffU2: [dxdx,dydx,dxdy,dydy]f
function DD=diffU2(x,y)
DDtol=1e-5;

    D1=(param(x+2*DDtol,y)-2*param(x+DDtol,y)+param(x,y))/(DDtol^2);
    D2=(param(x+DDtol,y+DDtol)-param(x,y+DDtol)-param(x+DDtol,y)+param(x,y))/(DDtol^2);
    D3=(param(x+DDtol,y+DDtol)-param(x+DDtol,y)-param(x,y+DDtol)+param(x,y))/(DDtol^2);
    D4=(param(x,y+2*DDtol)-2*param(x,y+DDtol)+param(x,y))/(DDtol^2);
    DD=[D1,D2,D3,D4];
    
end

% diffU3: [dxdxdx,dxdydx,dxdxdy,dxdydy,dydxdx,dydydx,dydxdy,dydydy]f

function DDD=diffU3(x,y)
DDDtol=1e-4;
    D1=(param(x+3*DDDtol,y)-3*param(x+2*DDDtol,y)+3*param(x+DDDtol,y)-param(x,y))/(DDDtol^3);
    D2=(param(x+2*DDDtol,y+DDDtol)-2*param(x+DDDtol,y+DDDtol)-param(x+2*DDDtol,y)+2*param(x+DDDtol,y)+param(x,y+DDDtol)-param(x,y))/(DDDtol^3);
    D3=(param(x+2*DDDtol,y+DDDtol)-param(x+2*DDDtol,y)-2*param(x+DDDtol,y+DDDtol)+2*param(x+DDDtol,y)+param(x,y+DDDtol)-param(x,y))/(DDDtol^3);
    D4=(param(x+DDDtol,y+2*DDDtol)-2*param(x+DDDtol,y+DDDtol)+param(x+DDDtol,y)-param(x,y+2*DDDtol)+2*param(x,y+DDDtol)-param(x,y))/(DDDtol^3);
    D5=(param(x+2*DDDtol,y+DDDtol)-2*param(x+DDDtol,y+DDDtol)+param(x,y+DDDtol)-param(x+2*DDDtol,y)+2*param(x+DDDtol,y)-param(x,y))/(DDDtol^3);
    D6=(param(x+DDDtol,y+2*DDDtol)-param(x,y+2*DDDtol)-2*param(x+DDDtol,y+DDDtol)+2*param(x,y+DDDtol)+param(x+DDDtol,y)-param(x,y))/(DDDtol^3);
    D7=(param(x+DDDtol,y+2*DDDtol)-2*param(x+DDDtol,y+DDDtol)-param(x,y+2*DDDtol)+2*param(x,y+DDDtol)+param(x+DDDtol,y)-param(x,y))/(DDDtol^3);
    D8=(param(x,y+3*DDDtol)-3*param(x,y+2*DDDtol)+3*param(x,y+DDDtol)-param(x,y))/(DDDtol^3);
    DDD=[D1,D2,D3,D4,D5,D6,D7,D8];
end

function MetricT=metricT(x,y)
    MetricT=zeros(2);
    D=diffU(x,y);
    for i=1:2
        for j=1:2
            MetricT(i,j)=D(:,i).'*D(:,j);
        end
    end
end

function DMetricTinv=dmetricTinv(x,y)
    G=metricT(x,y);
    dG=dmetricT(x,y);
    dG1=zeros(2);
    dG2=zeros(2);
    for i=1:2
        for j=1:2
            dG1(i,j)=dG(i,j);
            dG2(i,j)=dG(i,j+2);
        end
    end
    D1=-G^-1*dG1*G^-1;
    D2=-G^-1*dG2*G^-1;
    DMetricTinv=[D1,D2];
end

function DMetricT=dmetricT(x,y)
    DMetricT=zeros(2,4);
    D1=diffU(x,y);
    D2=diffU2(x,y);
    for k=0:1
        for i=1:2
            for j=1:2
                DMetricT(i,j+2*k)=D1(:,i).'*D2(:,j+2*k)+D2(:,i+2*k).'*D1(:,j);
            end
        end
    end
end

% ddmetricT: above: diff (2x4 mat) of dmetricT wrt x (l=1 upper, l=2 lower...), below diff of
% dmetricT wrt y
function DDMetricT=ddmetricT(x,y)
    DDMetricT=zeros(4,4);
    D1=diffU(x,y);
    D2=diffU2(x,y);
    D3=diffU3(x,y);
    for l=1:2
        for k=1:2
            for i=1:2
                for j=1:2
                    DDMetricT(i+2*(l-1),j+2*(k-1))=D3(:,i+2*(k-1)+4*(l-1)).'*D1(:,j)+D2(:,i+2*(k-1)).'*D2(:,j+2*(l-1))+D2(:,i+2*(l-1)).'*D2(:,j+2*(k-1))+D1(:,i).'*D3(:,j+2*(k-1)+4*(l-1));
                end
            end
        end
    end
end

function DChristoffel=dchristoffel(x,y)
    G=metricT(x,y);
    dGi=dmetricTinv(x,y);
    dG=dmetricT(x,y);
    ddG=ddmetricT(x,y);
    Gi=G^-1;
    DChristoffel=zeros(4);
    for l=1:2
        for k=1:2
            for i=1:2
                for j=1:2
                    DChristoffel(i+2*(l-1),j+2*(k-1))=1/2*(dGi(k,1+2*(l-1))*(dG(i,1+2*(j-1))+dG(j,1+2*(i-1))-dG(i,j))+dGi(k,2+2*(l-1))*(dG(i,2+2*(j-1))+dG(j,2+2*(i-1))-dG(i,j+2)))+1/2*(Gi(k,1)*(ddG(i+2*(l-1),1+2*(j-1))+ddG(j+2*(l-1),1+2*(i-1))-ddG(i+2*(l-1),j))+Gi(k,2)*(ddG(i+2*(l-1),2+2*(j-1))+ddG(j+2*(l-1),2+2*(i-1))-ddG(i+2*(l-1),j+2)));
                end
            end
        end
    end
end


%convention: k=1,2 is the first or second matrix (upper index in
%christoffel symbol, i,j are the usual matrix indices
function Christoffel=christoffel(x,y)
    G=metricT(x,y);
    Gi=G^-1;
    dG=dmetricT(x,y);
    Christoffel=zeros(2,4);
    for k=1:2
        for i=1:2
            for j=1:2
                Christoffel(i,j+2*(k-1))=1/2*(Gi(k,1)*(dG(i,1+2*(j-1))+dG(j,1+2*(i-1))-dG(i,j))+Gi(k,2)*(dG(i,2+2*(j-1))+dG(j,2+2*(i-1))-dG(i,j+2)));
            end
        end
    end
end

% riemannian curvature tensor: definiton: l in upper index, lower: ijk
function Rcurtensor=RcurT(x,y)
    Rcurtensor=zeros(4);
    dCS=dchristoffel(x,y);
    CS=christoffel(x,y);
    for l=1:2
        for k=1:2
            for i=1:2
                for j=1:2
                    Rcurtensor(i+2*(l-1),j+2*(k-1))=dCS(k+2*(i-1),j+2*(l-1))-dCS(k+2*(j-1),i+2*(l-1))+(CS(1,i+2*(l-1))*CS(k,j)-CS(1,j+2*(l-1))*CS(k,i)+CS(2,i+2*(l-1))*CS(k,j+2)-CS(2,j+2*(l-1))*CS(k,i+2));
                end
            end
        end
    end
end

function GaussCur=GC(x,y)
    Gi=metricT(x,y)^-1;
    R=RcurT(x,y);
    GaussCur=1/2*(Gi(1,1)*(R(1,1)+R(4,1))+Gi(2,1)*(R(1,2)+R(4,2))+Gi(1,2)*(R(1,3)+R(4,3))+Gi(2,2)*(R(1,4)+R(4,4)));
end

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
