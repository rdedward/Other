m=[1000]; % kg
J=[160]; % kgm^2
Zeta=[.5];
b1= 60; %50  60
Wn1=b1./(2.*m.*Zeta);
k1=((Wn1.^2).*m);
b2= 200; %70  200
Wn2=b2./(2.*m.*Zeta);
k2=((Wn2.^2).*m); 
w1=12.8803; %radians/seconds
w2=20.4202; %radians/seconds
l1=1; % m
l2=1; % m
t=0;
top1=(m.*((w1).^2)+k1+k2-b1.*(w1).*tan(w1.*t+atan(-2.*(.5).*(w1./sqrt(k1./m)./(1-((w1./sqrt(w1./sqrt(k1./m)))))))));
bot1=(-k1.*l1+k2.*l2+b1.*l1.*w1.*tan(w1.*t+atan(-2.*(.5).*(w1./sqrt(k1./m)./(1-((w1./sqrt(w1./sqrt(k1./m)))))))));
top2=(m.*((w2).^2)+k1+k2-b1.*(w2).*tan(w2.*t+atan(-2.*(.5).*(w2./sqrt(k1./m)./(1-((w2./sqrt(w2./sqrt(k1./m)))))))));
bot2=(-k1.*l1+k2.*l2+b1.*l1.*w2.*tan(w2.*t+atan(-2.*(.5).*(w2./sqrt(k1./m)./(1-((w2./sqrt(w2./sqrt(k1./m)))))))));
r1=top1./bot1
r2=top2./bot2


%% Part B


K1 = k1;
K2 = k2;
B1 = b1;
B2 = b2;
w = sym('w');
t = sym('t');
m = m;
phi = atan(-2.*(.5).*(w./sqrt(k1./m)./(1-((w./sqrt(w./sqrt(k1./m)))))));



L1 = l1;
L2 = l2;
J = J;
M = [-m*w^2+K1+K2-B1*w*tan(w*t+phi)-B2*w*tan(w*t+phi),-K1*L1+K2*L2+B1*L1*w*tan(w*t*phi)-B2*L2*w*tan(w*t+phi);-K1*L1+K2*L2+B1*L1*w*tan(w*t*phi)-B2*L2*w*tan(w*t+phi),-J*w^2+K1*L1^2+K2*L2^2-B1*L1^2*w*tan(w*t+phi)+B2*L2^2*w*tan(w*t+phi)]

D = det(M)

yf = .15.*sin(733.*t);
yfdot = 110.*cos(733.*t);
yr = .15.*sin(733.*t-2.51);
yrdot = 110.*cos(733.*t-2.51);

Fo1 = -k1.*yr-k2.*yf-b1.*yrdot-b2.*yfdot;
Fo2 = (k1*l1*yr-k2.*l2*yf)+(b1.*l1.*yrdot-b2.*l2.*yfdot);

Xbar = (M(1,1).*Fo1+M(1,2).*Fo2)/D;
Thetabar = (M(2,1).*Fo1+M(2,2).*Fo2)/D;

%For First frequency
Xbarsub = subs(Xbar,w,w1)
Thetabarsub = subs(Thetabar,w,w1)

ezplot(Xbarsub,[0,1])
figure;
ezplot(Thetabarsub,[0,1])

%For Second frequency
Xbarsub2 = subs(Xbar,w,w2)
Thetabarsub2 = subs(Thetabar,w,w2)
figure;
ezplot(Xbarsub2,[0,1])
figure;
ezplot(Thetabarsub2,[0,1])







