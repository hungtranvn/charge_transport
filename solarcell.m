function Ia = solarcell(Va, G, TaC, Rs, Rsh)

    k = 1.381e-23; % Boltzmann constant
    q = 1.602e-19;% Electron charge

    n = 1.7; % Diode ideality factor (n),
    Eg = 2.3; % Band gap energy 1.8eV PCBM:P3HT (2.6?)
    TrK = 298; % Reference temperature (25C) in Kelvin
    Voc_TrK = 0.58; % Voc (open circuit voltage per cell) @ temp TrK
    Isc_TrK = 6.12e-3; % Isc (short circuit current per cell) @ temp TrK
    a = 5.1e-3; % Temperature coefficient of Isc (0.065%/C)
    TaK = 273 + TaC; % Module temperature in Kelvin
    Vc = Va; % Cell voltage

    % Calculate short-circuit current for TaK
    Isc = Isc_TrK * (1 + (a * (TaK - TrK)));

    % Calculate photon generated current @ given irradiance
    Iph = G * Isc;

    % Define thermal potential (Vt) at temp TrK
    Vt_TrK = n * k * TrK / q;

    % Define b = Eg * q/(n*k);
    b = Eg * q /(n * k);

    % Calculate reverse saturation current for given temperature
    Ir_TrK = Isc_TrK / (exp(Voc_TrK / Vt_TrK) -1);
    Ir = Ir_TrK * (TaK / TrK)^(3/n) * exp(-b * ((1 / TaK) - (1 / TrK)));

    % Define thermal potential (Vt) at temp Ta
    Vt_Ta = n * k * TaK / q;

    Ia = zeros(size(Vc )); % Initialize Ia with zeros

    % Solve for Ia by Newton's method: Ia2 = Ia1 - f(Ia1)/f'(Ia1)
    % Perform 5 iterations
    for j=1:10
        Ia = Ia  + (-Iph - Ia + Ir .* ( exp((Vc - Ia .* Rs) ./ Vt_Ta) -1) + (Vc - Ia.*Rs)/(Rsh))./ (1 + Ir * (Rs / Vt_Ta) .* exp((Vc - Ia .* Rs) ./ Vt_Ta) + Rs/Rsh);
end
