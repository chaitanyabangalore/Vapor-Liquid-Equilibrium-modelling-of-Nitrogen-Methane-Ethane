

function binary
global Tc R b ac kappa
Tc = [305.4 190.4];
Pc = [48.8 46];
R = 83.14; 
w = [0.099 0.011];
options=optimset('Display','off');


b = 0.0778*R*Tc./Pc;
ac = 0.45724*R^2*Tc.^2./Pc;
kappa = 0.37464+1.54226*w-0.269932*w.^2;
T=352; %K

pm=0.055;
y=0.6;
xx=0:0.2:1;
nx=length(xx);
for i=1:nx
    x=xx(i);

    Vl=1.5*b(1);
    Vv=R*T/pm;
    pold=100;

    while abs(pold-pm)>0.00001
        Vl=fsolve(@(Vl) pm-p(T,Vl,x),Vl,options);
        phil1=exp(lnphi1(pm,Vl,T,x));
        phil2=exp(lnphi2(pm,Vl,T,x));
        yold=2;

        while abs(y-yold)>0.00001
            Vv=fsolve(@(Vv) pm-p(T,Vv,y),Vv,options);
            phiv1=exp(lnphi1(pm,Vv,T,y));
            phiv2=exp(lnphi2(pm,Vv,T,y));

            K1=phil1/phiv1;
            K2=phil2/phiv2;

            sum1=K1*x+K2*(1-x);
            yold=y;
            y=K1*x;

        end
   

    pold=pm;
    pm=pm*sum1;

    end
    pp(i)=pm;
    yy(i)=y;
end

pp;
yy;
hold on
plot(xx,pp, 'b')
plot(yy,pp, 'r')
end

function y=am(x,t)
global Tc kappa ac b
alpha=(1+kappa.*(1-sqrt(t./Tc))).^2;
a=ac.*alpha;
bmix = x*b(1)+(1-x)*b(2);
kij=-0.0026;
y=x^2.*a(1)+(1-x)^2.*a(2)+2*x*(1-x)*sqrt(a(1).*a(2))*(1-kij);
end

function y=p(t,v,x)
global R b
bmix = x*b(1)+(1-x)*b(2);
y=R*t./(v-bmix)-am(x,t)./(v*(v+bmix)+bmix*((v-bmix)));
end

function y=lnphi1(p,v,t,x)
global R b ac kappa Tc
s2 = sqrt(2);
z=p*v/R/t;
alpha=(1+kappa.*(1-sqrt(t./Tc))).^2;
a=ac.*alpha;
bmix = x*b(1)+(1-x)*b(2);
y = b(1)/bmix*(z-1)-log(z-bmix*p/R/t)-am(x,t)/2/s2/bmix/R/t*(2*(x*a(1)+(1-x)*sqrt(a(1)*a(2)))/am(x,t)-b(1)/bmix)*log((z+(1+s2)*bmix*p/R/t)/(z+(1-s2)*bmix*p/R/t));
end

function y=lnphi2(p,v,t,x)
global R b ac kappa Tc 
s2 = sqrt(2);
z=p*v/R/t;
alpha=(1+kappa.*(1-sqrt(t./Tc))).^2;
a=ac.*alpha;
bmix = x*b(1)+(1-x)*b(2);
y=b(2)/bmix*(z-1)-log(z-bmix*p/R/t)-am(x,t)/2/s2/bmix/R/t*(2*(x*sqrt(a(1)*a(2))+(1-x)*a(2))/am(x,t)-b(2)/bmix)*log((z+(1+s2)*bmix*p/R/t)/(z+(1-s2)*bmix*p/R/t));
end

