%Joseph Filbert
%EE515
%Project

cc=299792458; %m/s
mu0 = 4*pi*1e-7; %H/m
eps0 = 1/cc^2/mu0; %F/m

%waveguide characteristic impedance Z_TE
a = 0.900*25.4e-3; %m
b = 0.400*25.4e-3; %m
Zte = @(f) 2*pi*mu0*f./sqrt((2*pi*f/cc).^2-(pi/a)^2);
w = linspace(a/2,a/10,10);
for ii=1:10
%load S parameters
fName = sprintf('EE515waveguideFilterProject_%d_%d.s2p',ii,ii);
sobj=sparameters(fName);
gamma = squeeze(sobj.Parameters(1,1,:));
Zin = Zte(sobj.Frequencies).*(1+gamma)./(1-gamma);

if ii==1
fig1 = figure; hold on; fig1.Color='white';
ax1=gca;
fig2 = figure; hold on; fig2.Color='white';
ax2=gca;
end
p1 = plot(ax1,sobj.Frequencies./1e9,real(Zin),'DisplayName',sprintf('%.4f',w(ii)));
p2 = plot(ax1,sobj.Frequencies./1e9,imag(Zin)./(2*pi*sobj.Frequencies)./1e-9);
p2.HandleVisibility='off';
p2.Color = p1.Color;
p2.LineStyle = '-.';

powerLoss = 10*log10(squeeze(1-abs(sobj.Parameters(1,1,:)).^2+abs(sobj.Parameters(2,1,:)).^2));
p3 = plot(ax2,sobj.Frequencies/1e9,powerLoss,'DisplayName',sprintf('%.4f',w(ii)));
p3.Color=p1.Color;
end
ax1.XLabel.String='Frequency (GHz)';
ax1.YLabel.String= 'Z (\Omega), L (nH)';
ldg=legend(ax1);ldg.Title.String='w (mm)';
ldg.Location='north';ldg.NumColumns=5;

ax2=gca;ax2.XLabel.String='Frequency (GHz)';
ax2.YLabel.String='(1-|S_{11}|^2-|S_{21}|^2) dB';
ldg = legend(ax2);ldg.Title.String='w (mm)';
ldg.Location='south';ldg.NumColumns=5;
%this works for a single iris width, need to get the other iris widths, and
%then do simulation for 