clear all

%givens
mdot = 50.33; %kg/s
cz = 288; %m/s
a1 = -15; %deg
To1 = 1910; %K
Po1 = 1.293e6; %Pa
MW = 28.8; %g/mol
gamma = 1.29;
mu = 6.52e-5; %Ns/m^2
R = (8314/MW); %J/kgK
Cp = (gamma*R)/(gamma-1); %J/kgK
etaTR = 0.99;
Wproportion = 0.675;
N = 7920; %rpm

%constraints and requirements
M2 = 1;
a3min = -45; %deg
a3max = 45; %deg
zweif = 0.8;
Remin = 1e5;
Remax = 5e5;
Umax = 423; %m/s - specifically at blade tip
Wout = 20.6e6 * Wproportion; %W

%lower and upper bounds
lb = [a3min Remin Remin];
ub = [a3max Remax Remax];

%starting points
a3guess = -6; %deg - this should probably go to -11 degrees to max
Reonguess = 1.49e5; %reynolds number at nozzle throat
Reorguess = 1.13e5; %reynolds number at rotor "throat"
x0 = [a3guess Reonguess Reorguess];

%perform optimization
objectiveFunction = @(x) optimizeStage(x, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout);
nonlcon = @(x) constraints(x, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout);
[out, etaST, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon);
etaST = -etaST; %flip to positive since we're done tricking fmincon

%postprocess and output the necessary data for the homework
postProcess(out, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout);

%stage efficiency optimization function
function etaST = optimizeStage(x, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout)
    %unpack design variables
    [a3, Reon, Reor] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, c2, ctheta2, a2, w2, wtheta2, b2, w3, wtheta3, b3, c3, ctheta3, M1] = velocityTriangleAnalysis(Wout, N, R, cz, a1, a3, Po1, To1, M2, Cp, gamma, mdot, etaTR);

    %calculate nozzle radii and pitchline U
    [rtn, rhn, rm, U, rho1] = nozzleRadiiAndU(Wout, N, gamma, R, Po1, cz, a1, To1, Cp, mdot, ctheta2, ctheta3, etaTR);
    
    %calculate psi, phi, Rdeg
    [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3);

    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigmarz, sigman, sigmar] = sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rtn - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rm, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1);

    %calculate blade aspect ratio and pressure ratio for rotor
    To2 = To1;
    [rtr, rhr] = rotorRadii(cz, rm, mdot, To2, Po2, M2, gamma, R);

    h2 = rtr - rhr;
    [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar] = rotor_blade_aspect_ratio(U, sigmar, sigmarz, rm, Reor, mu, R, Cp, To1, Po1, Po2, gamma, M2, c2, w2, w3, h2, c3, psi, b2, b3);

    num = psi*(gamma-1)*((U^2)/(gamma*R*To1));
    den = Prtot^((1-gamma)/gamma) - 1;
    etaST = -num/den; %negative so fmincon will maximize
end

%nonlinear constraint function
function [c, ceq] = constraints(x, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout)
    %unpack design variables
    [a3, Reon, Reor] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, c2, ctheta2, a2, w2, wtheta2, b2, w3, wtheta3, b3, c3, ctheta3, M1] = velocityTriangleAnalysis(Wout, N, R, cz, a1, a3, Po1, To1, M2, Cp, gamma, mdot, etaTR);

    %calculate nozzle radii and pitchline U
    [rtn, rhn, rm, U, rho1] = nozzleRadiiAndU(Wout, N, gamma, R, Po1, cz, a1, To1, Cp, mdot, ctheta2, ctheta3, etaTR);

    %calculate psi, phi, Rdeg
    [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3);

    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigmarz, sigman, sigmar] = sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rtn - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rm, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1);

    %calculate rotor radii
    To2 = To1;
    [rtr, rhr] = rotorRadii(cz, rm, mdot, To2, Po2, M2, gamma, R);

    %constrain blade tip speed
    c(1) = ((rtr * pi * N)/30) - Umax;
    ceq = [];
end

%work calculation function
function Wdot = calculateWorkDone(etaTR, mdot, U, ctheta3, ctheta2)
    Wdot = -(etaTR * mdot * U * (ctheta3 - ctheta2)); %multiplying by negative one, making work done > 0
end

%unpacking function
function [a3, Reon, Reor] = unpack(x)
    a3 = x(1);
    Reon = x(2);
    Reor = x(3);
end

%solidity calculation function
function [sigmanz, sigmarz, sigman, sigmar] = sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma)
    %solve for nozzle solidity
    anmean = atand(((ctheta1+ctheta2)/2)/cz); % stagger angle
    sigmanz = ((ctheta1/ctheta2 - 1)*sind(2*a2)*(((gamma/2)*(M2^2))/((1 + ((gamma-1)/2)*(M2^2))^(gamma/(gamma-1)) - 1)))/(-zweif); % axial solidity
    sigman = sigmanz/cosd(anmean); % regular solidity

    %calculate M3rel
    To2 = To1;
    To3 = To2 + (psi * U^2)/Cp;
    T3 = To3 - (c3^2)/(2*Cp);
    M3 = c3/sqrt(gamma*R*T3);
    M3rel = (M3/c3)*w3;
     
    %solve for rotor solidity
    brmean = atand(((wtheta2+wtheta3)/2)/cz); % stagger angle
    sigmarz = ((wtheta2/wtheta3 - 1)*sind(2*b3)*(((gamma/2)*(M3rel^2))/((1 + ((gamma-1)/2)*(M3rel^2))^(gamma/(gamma-1)) - 1)))/(zweif); % axial solidity
    sigmar = sigmarz/cosd(brmean); % regular solidity
end

%coefficient calculation function
function [rtn, rhn, rm, U, rho1] = nozzleRadiiAndU(Wout, N, gamma, R, Po1, cz, a1, To1, Cp, mdot, ctheta2, ctheta3, etaTR)
    %calculate M1
    c1 = cz/cosd(a1);
    T1 = To1 - (c1^2)/(2*Cp);
    M1 = c1/sqrt(gamma*R*T1);

    %calculate necessary rm using wdot requirement
    rm = (30*Wout/(pi*etaTR*mdot*N*(ctheta2-ctheta3)));
    P1 = Po1 * (1 + ((gamma-1)/2)*(M1^2))^(-gamma/(gamma-1));
    T1 = To1 * (1 + ((gamma-1)/2)*(M1^2))^(-1);
    rho1 = P1/(R*T1);
    rhn = sqrt(rm^2 - (mdot/(2*rho1*cz*pi))); % hub radius
    rtn = sqrt(2*rm^2 - rhn^2); % tip radius
    U = (rm*N*pi/30); % blade velocity at mean radius    
end

function [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3)
    phi = cz/U; % flow coefficient phi
    psi = (ctheta3-ctheta2)/U; % stage loading coefficient psi
    Rdeg = (-phi/2)*(tand(b2) + tand(b3)); % degree of reaction
end

%rotor tip and hub radius calculation function
function [rtr, rhr] = rotorRadii(cz, rm, mdot, To2, Po2, M2, gamma, R)
    %calculate flow properties at rotor inlet
    P2 = Po2 * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
    T2 = To2 * (1 + ((gamma-1)/2)*(M2^2))^(-1);
    rho2 = P2/(R*T2);

    %apply continuity and assume constant pitchline radius to obtain radii
    rtr = sqrt((mdot/(2*rho2*cz*pi)) + rm^2);
    rhr = sqrt(2*(rm^2) - rtr^2);
end

%velocity triangle analysis function
function [c1, ctheta1, c2, ctheta2, a2, w2, wtheta2, b2, w3, wtheta3, b3, c3, ctheta3, M1] = velocityTriangleAnalysis(Wout, N, R, cz, a1, a3, Po1, To1, M2, Cp, gamma, mdot, etaTR)
    %nozzle inlet
    c1 = cz/cosd(a1);
    ctheta1 = cz*tand(a1);
    T1 = To1 - (c1^2)/(2*Cp);
    M1 = c1/sqrt(gamma*R*T1);

    %nozzle outlet
    To2 = To1; %assuming adiabatic across nozzle
    T2 = To2*(1 + ((gamma-1)/2)*M2^2)^-1;
    c2 = M2 * sqrt(gamma*R*T2);
    a2 = acosd(cz/c2);
    ctheta2 = cz*tand(a2);

    %rotor outlet - stationary frame
    ctheta3 = cz*tand(a3);
    c3 = cz/cosd(a3);

    %determine U from work requirement
    [rtn, rhn, rm, U, rho1] = nozzleRadiiAndU(Wout, N, gamma, R, Po1, cz, a1, To1, Cp, mdot, ctheta2, ctheta3, etaTR);

    %rotor outlet - rotating frame
    wtheta3 = ctheta3 - U;
    w3 = sqrt(wtheta3^2 + cz^2);
    b3 = atand(wtheta3/cz);

    %rotor inlet
    wtheta2 = ctheta2 - U;
    w2 = sqrt(wtheta2^2 + cz^2);
    b2 = atand(wtheta2/cz);

end

function [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar] = rotor_blade_aspect_ratio(U, sigmar, sigmarz, rm, Reor, mu, R, Cp, To1, Po1, Po2, gamma, M2, c2, w2, w3, h2, c3, psi, b2, b3)
    iterations = 2000; % iterations for both loops
    To2 = To1;
    T2 = To2 - (c2^2)/(2*Cp);
    To3 = To2 + (psi * U^2)/Cp;
    %transform into rotor reference frame
    M2rel = (M2/c2)*w2;
    Po2rel = Po2*((1+((gamma-1)/2)*M2rel^2)^(gamma/(gamma-1)))/((1+((gamma-1)/2)*M2^2)^(gamma/(gamma-1)));
    To2rel = T2*(1+((gamma-1)/2)*M2rel^2);
    PoIterable = Po2rel;

    T3 = To3 - (c3^2)/(2*Cp);
    M3 = c3/sqrt(gamma*R*T3);
    M3rel = (M3/c3)*w3;

    %blade aspect ratio - rotor
    for i = 1:iterations
        T3 = To3 * (1 + ((gamma-1)/2)*(M3^2))^-1;
        P3 = PoIterable(i) * (1 + ((gamma-1)/2)*(M3rel^2))^(-gamma/(gamma-1));
        rho3 = P3/(R*T3);
        mur = mu;
        nur = mur/rho3;
        
        or = Reor * nur / w3;
        sr = or / cosd(b3);
        Nr = getIntegerBladeCount(2*pi*rm/sr); %num of blades for rotor
        sr = (2*pi*rm)/Nr;
        br = sigmar*sr;
        bladeARr = h2/br;
        
        %stagnation pressure ratio - rotor
        bzr = sigmarz * sr;
        Crotor = 0.975 + 0.075*(bzr/h2);
        
        zetastarr = 1.04 + 0.06*((b2 + b3)/100)^2;
        Dhr = ((2*sr*h2*cosd(b3))/(sr*cosd(b3)+h2));
        Rer = (rho3 * w3 * Dhr) / mur;
        zetar = (zetastarr*Crotor - 1)*((10^5)/Rer)^0.25;
        
        num = (1 - ((w3^2)/(2*Cp*To2rel)) * (1/(1-zetar)));
        den = (1 - ((w3^2)/(2*Cp*To2rel)));
        Po3rel = Po2rel - ((1 - (num/den)^(gamma/(gamma-1)))*Po2rel);
        
        PoIterable(i+1) = Po3rel;
    end
    %convert back into absolute reference frame
    Po3 = Po3rel*((1+((gamma-1)/2)*M3^2)^(gamma/(gamma-1)))/((1+((gamma-1)/2)*M3rel^2)^(gamma/(gamma-1)));
    Prtot = Po1/Po3;
end
 
function [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rm, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1)
    iterations = 2000; % iterations for both loops
    %blade aspect ratio - nozzle
    PoIterable = Po1;
    for i = 1:iterations
        T2 = To1 * (1 + ((gamma-1)/2)*(M2^2))^-1;
        P2 = PoIterable(i) * (1 + ((gamma-1)/2)*(M2^2))^(-gamma/(gamma-1));
        rho2 = P2/(R*T2);
        
        mun = mu;
        nun = mun/rho2;
        
        c2 = sqrt(ctheta2^2 + cz^2);
        on = Reon * nun / c2;
        sn = on / cosd(a2);
        Nn = getIntegerBladeCount(2*pi*rm/sn); %num of blades for nozzle
        sn = (2*pi*rm)/Nn;
        bn = sigman*sn;
        bladeARn = h1/bn;
        
        %stagnation pressure ratio - nozzle
        bzn = sigmanz * sn;
        Cnozzle = 0.993 + 0.021*(bzn/h1);
        
        zetastarn = 1.04 + 0.06*((a1 + a2)/100)^2;
        Dhn = ((2*sn*h1*cosd(a2))/(sn*cosd(a2)+h1));
        Ren = (rho2 * c2 * Dhn) / mun;
        zetan = (zetastarn*Cnozzle - 1)*((10^5)/Ren)^0.25;
        
        num = (1 - ((c2^2)/(2*Cp*To1)) * (1/(1-zetan)));
        den = (1 - ((c2^2)/(2*Cp*To1)));
        Po2 = Po1 - ((1 - (num/den)^(gamma/(gamma-1)))*Po1);
        PoIterable(i+1) = Po2;
    end
end
 
%blade count rounding function
function Nrounded = getIntegerBladeCount(Nunrounded)
    Nrounded = floor(Nunrounded);
    if (mod(Nrounded,2) == 0)
        Nrounded = Nrounded + 1;
    end
end

%post-processing and output function
function postProcess(x, R, To1, Po1, M2, cz, a1, Cp, gamma, Umax, N, mdot, zweif, mu, etaTR, Wout)
    %unpack design variables
    [a3, Reon, Reor] = unpack(x);

    %perform velocity triangle analysis with chosen parameters
    [c1, ctheta1, c2, ctheta2, a2, w2, wtheta2, b2, w3, wtheta3, b3, c3, ctheta3, M1] = velocityTriangleAnalysis(Wout, N, R, cz, a1, a3, Po1, To1, M2, Cp, gamma, mdot, etaTR);

    %calculate nozzle radii and pitchline U
    [rtn, rhn, rm, U, rho1] = nozzleRadiiAndU(Wout, N, gamma, R, Po1, cz, a1, To1, Cp, mdot, ctheta2, ctheta3, etaTR);
    
    %calculate psi, phi, Rdeg
    [psi, phi, Rdeg] = psiPhiR(cz, U, ctheta3, ctheta2, b2, b3);

    %calculate optimum solidities using zweifel criteria
    [sigmanz, sigmarz, sigman, sigmar] = sigmaZweif(w3, To1, U, Cp, R, cz, ctheta1, ctheta2, wtheta2, wtheta3, a2, b3, M2, psi, c3, zweif, gamma);

    %calculate blade aspect ratio and pressure ratio for nozzle
    h1 = rtn - rhn;
    [bladeARn, Po2, Nn, zetan] = nozzle_blade_aspect_ratio(rm, Reon, mu, R, Cp, To1, Po1, gamma, M2, sigman, sigmanz, a1, a2, ctheta2, cz, h1);

    %calculate blade aspect ratio and pressure ratio for rotor
    To2 = To1;
    [rtr, rhr] = rotorRadii(cz, rm, mdot, To2, Po2, M2, gamma, R);

    h2 = rtr - rhr;
    [bladeARr, Prtot, To2, To3, rho3, Nr, M2rel, zetar] = rotor_blade_aspect_ratio(U, sigmar, sigmarz, rm, Reor, mu, R, Cp, To1, Po1, Po2, gamma, M2, c2, w2, w3, h2, c3, psi, b2, b3);

    num = psi*(gamma-1)*((U^2)/(gamma*R*To1));
    den = Prtot^((1-gamma)/gamma) - 1;
    etaST = num/den;

    %perform final calculations needed and prepare data for output

    %for general outputs
    Po3 = Po1/Prtot;
    num = log(To1/To3);
    den = log(1 + (1/etaST)*(To1/To3 - 1));
    etaPoly = num/den;
    rho3_rho1 = rho3/rho1;

    row1 = {'', '°R', 'U (m/s)', 'rm (m)', 'phi', 'psi', 'To3 (K)', 'Po3 (kPa)', 'etaST (%)', 'etaPoly (%)', 'rho3/rho1', '', ''};
    row2 = {'', Rdeg, U, rm, phi, psi, To3, (Po3/1000), (etaST*100), (etaPoly*100), rho3_rho1, '', ''};

    %for nozzle outputs
    chi2 = a2 + 0; %assuming no deviation since M2=1
    Po2_Po1 = Po2/Po1;
    
    row3 = {'', 'chi2 (deg)', 'a2 (deg)', 'sigmaN', 'Nn', 'Reon', 'rhn', 'rtn', 'bladeARn', 'zetan', 'Po2/Po1', 'M1', ''};
    row4 = {'', chi2, a2, sigman, Nn, Reon, rhn, rtn, bladeARn, zetan, Po2_Po1, M1, ''};

    %for rotor outputs
    chi3 = b3 - cartersDeviation((b3-b2), sigmar); %CHECK
    Po3_Po2 = (1/Prtot) * (1/Po2_Po1);
    AN2 = pi*(rtr^2 - rhr^2) * (N^2);

    row5 = {'b3 (deg)', 'chi3 (deg)', 'a3 (deg)', 'sigmaR', 'Nr', 'Reor', 'rhr', 'rtr', 'bladeARr', 'zetar', 'Po3/Po2', 'M2rel', 'AN^2 (m^2 RPM^2)'};
    row6 = {b3, chi3, a3, sigmar, Nr, Reor, rhr, rtr, bladeARr, zetar, Po3_Po2, M2rel, AN2};

    %concatenate output matrix
    output = [row1; row2; row3; row4; row5; row6];
    
    %write file
    filename = 'turbineData.xlsx';
    writecell(output, filename);
end

function d = cartersDeviation(flowAngleDelta, sigma)
    m = 1/8;
    d = m*(abs(flowAngleDelta)/sigma);
end