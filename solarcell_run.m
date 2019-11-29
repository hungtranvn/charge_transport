% Define constant
G = 1;
TaC = 25;
Rs = 5;
Rsh = 900;

% Functions to plot
figure
hold on

Vc = [-0.1:0.02:0.7];
Ic = solarcell(Vc, G, TaC, Rs, Rsh);
plot(Vc, Ic,'marker','*');
title('Photovoltaic I-V Curve')
xlabel('Cell Voltage (V)')
ylabel('Cell Current (A)')
axis([-0.1 0.8 -0.012 0.001])
h = [Vc;Ic];
