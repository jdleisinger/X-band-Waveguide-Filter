%Joseph Filbert
cc = 299792458; mu0=4*pi*1e-7; eps0 = 1/cc^2/mu0;
a=0.900*25.4e-3;
f0 = 10e9;
BW = 200e6;
fcH = f0+BW/2;
fcL = f0-BW/2;
Ap = 3; %dB
N = 5;
Rp = [0.5, 3];

robj = rffilter();
robj=repmat(robj,1,3);
ZTE10 = @(w) (w/cc)*sqrt(mu0/eps0)/sqrt((w/cc)^2-(pi/a)^2);
%%
robj(1) = rffilter('ResponseType','Bandpass','Implementation','LC Pi','PassbandFrequency',[fcL fcH],                        ...
    'PassbandAttenuation',Ap,'FilterOrder',N,'Zin',ZTE10(2*pi*f0),'Zout',ZTE10(2*pi*f0));
robj(1).Name = 'Butterworth';
robj(2) = rffilter('FilterType','Chebyshev','ResponseType','Bandpass','Implementation','LC Pi','PassbandFrequency',[fcL fcH],                        ...
    'PassbandAttenuation',Rp(1),'FilterOrder',N,'Zin',ZTE10(2*pi*f0),'Zout',ZTE10(2*pi*f0));
robj(2).Name = 'Chebyshev0p5dB';
robj(3) = rffilter('FilterType','Chebyshev','ResponseType','Bandpass','Implementation','LC Pi','PassbandFrequency',[fcL fcH],                        ...
    'PassbandAttenuation',Rp(2),'FilterOrder',N,'Zin',ZTE10(2*pi*f0),'Zout',ZTE10(2*pi*f0));
robj(3).Name = 'Chebyshev3p0dB';
%%
freq = linspace(8.2,12.4,1001)*1e9;
figure; hold on; fig=gcf; fig.Color='white';
for ii=1:3
s=sparameters(robj(ii),freq,ZTE10(2*pi*f0));
line = rfplot(s,2,1); line.DisplayName = robj(ii).Name;
rfwrite(s,strcat('EE515_',robj(ii).Name))
end
lgd = legend;ldg.Title='S_{21} dB';
lgd.Location = "best";

%%
syms k l c w l1 l2
abcdC=[0, 1i*k; 1i/k, 0]*[1 , 1i*w*l; 0, 1]*[0, 1i*k; 1i/k, 0]*[-1,0;0,-1]
abcdE=[1, 1i*w*(l1-l); 0, 1]*[1, 0; 1i*w*c, 1]*[1, 1i*w*(l2-l/2); 0, 1]
abcdD=[0, 1i*k; 1i/k, 0]*[1 , 1i*(w^2*l*c-1)/w/c; 0, 1]*[0, 1i*k; 1i/k, 0]*[-1,0;0,-1]
%%
w1=2*pi*fcL; w2=2*pi*fcH;
w0 = 2*pi*f0;
deltaW = (w2-w1);
fracBW = deltaW/w0;
Zte = ZTE10(2*pi*f0);
g0 = 1.0000; g1 = 1.7058; g2 = 1.2296; 
g3 = 2.5408; g4 = g2; g5 = g1; g6 = 1.0000;
ff = pi/2;
% KS1 = Zte*sqrt(Zte*deltaW*ff/g0/g1);
% K12 = deltaW*ff*sqrt(1/g1/g2);
% K23 = deltaW*sqrt(1/g2/g3);
% K34 = deltaW*sqrt(1/g3/g4);
% K45 = deltaW*sqrt(1/g4/g5);
% K5L = sqrt(Zte*deltaW/g5/g6);

KS1 = Zte*sqrt(fracBW*ff/g0/g1);
K12 = Zte*fracBW*ff*sqrt(1/g1/g2);
K23 = Zte*fracBW*ff*sqrt(1/g2/g3);
K34 = Zte*fracBW*ff*sqrt(1/g3/g4);
K45 = Zte*fracBW*ff*sqrt(1/g4/g5);
K5L = Zte*sqrt(fracBW*ff/g5/g6);


XS1 = KS1/(1-(KS1/Zte)^2); thetaS1 = -atan(2*XS1/Zte);
X12 = K12/(1-(K12/Zte)^2); theta12 = -atan(2*X12/Zte);
X23 = K23/(1-(K23/Zte)^2); theta23 = -atan(2*X23/Zte);
X34 = K34/(1-(K34/Zte)^2); theta34 = -atan(2*X34/Zte);
X45 = K45/(1-(K45/Zte)^2); theta45 = -atan(2*X45/Zte);
X5L = K5L/(1-(K5L/Zte)^2); theta5L = -atan(2*X5L/Zte);

lambdaG = 2*pi/sqrt((w0/cc)^2-(pi/a)^2);

nXS1=XS1*lambdaG/Zte/2/a;
nX12=X12*lambdaG/Zte/2/a; L12 = lambdaG/2*(1+(thetaS1+theta12)/2/pi);
nX23=X23*lambdaG/Zte/2/a; L23 = lambdaG/2*(1+(theta12+theta23)/2/pi);
nX34=X34*lambdaG/Zte/2/a; L34 = lambdaG/2*(1+(theta23+theta34)/2/pi);
nX45=X45*lambdaG/Zte/2/a; L45 = lambdaG/2*(1+(theta34+theta45)/2/pi);
nX5L=X5L*lambdaG/Zte/2/a; L5L = lambdaG/2*(1+(theta45+theta5L)/2/pi);