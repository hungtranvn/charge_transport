function [ rvals ] = Ortizconde( IVData )
    %Constants
    q = 1.602*10^-19;
    k = 1.38*10^-23;
    Tc = 298;
    Vth = k*Tc/q;
    n = 1.5;
    Pin = 100; %mW/cm^2
    %Breaking apart input matrix
    V=IVData(:,1);
    I=-IVData(:,2);

    %Finding the Isc (foo is a dummy variable)
    [foo, SCindex] = min(abs(V));
    Isc = I(SCindex);

    %Finding the Voc
    [foo, OCindex] = min(abs(I));
    Voc = V(OCindex);

    V = IVData(SCindex:OCindex,1); 
    I = -IVData(SCindex:OCindex,2);
    Isc_I = Isc - I;
    
    y = quadraticsplinesintegral(V, Isc_I);
    V = IVData(SCindex+1:OCindex,1);
    I = -IVData(SCindex+1:OCindex,2);
    Isc_I = Isc - I;
    Isc_I = -Isc_I;
    cof = multiregression(V, Isc_I, y);
    
    %Finding the FF(fill factor)
    V=IVData(:,1);
    I=IVData(:,2);
    Vmod = V(SCindex:OCindex);
    Imod = I(SCindex:OCindex);

    P = Vmod.*Imod;
    [Pmax,Pmaxpt] = max(abs(P));
    Vprime = Vmod(Pmaxpt);
    Iprime = Imod(Pmaxpt);

    FF = -(Iprime.*Vprime)./(Voc.*Isc);
    
    %Outputting relevant data
    
    PCE = FF.*Isc.*Voc./0.106;
    % Estimate Rsh, Rs, I0, Il values
    Gp = 2*cof(4);
    Rsh = 1/Gp;
    n1 = cof(1)*((sqrt(1 + 16*cof(4)*cof(5)) - 1) + 4*cof(2)*cof(4))/(4*Vth*cof(4));
    Rs = (sqrt(1+16*cof(4)*cof(5)) - 1)/(4*cof(4));
    I0 = (Isc - Voc/Rsh )*exp(-Voc/(n*Vth));
    Il = Isc - I0 + I0*exp((Rs*Isc)/(n*Vth)) + (Rs*Isc)/Rsh;
    
    rvals = [Rs; Rsh; FF; PCE];
