%params
Vabs = 2;
L = 200e-6;
Rga = 0.1;
inj = 0.6;


Pin  = 1e-3 ;%10e-3; %min Pin with I_bias 24.45e-3 and pulse_duration 0.9e-9
Pulse_duration = 0.9e-9;

I_bias = 24.35;

generatePhasePortrait();


function generatePhasePortrait()
    % Define parameters for the laser model
    Vabs = 1;  % Example value for Vabs (adjust as needed)
    L = 1e-3;  % Example value for L (adjust as needed)
    Rga = 0.1; % Example value for Rga (adjust as needed)
    Pin = 1e-3; % Example value for Pin (adjust as needed)
    Pulse_duration = 1e-9; % Example value for Pulse_duration (adjust as needed)

    % Get constants and initial conditions
    p = constants(Vabs, L, Rga);
    [Ng, Nabs, Nph, Pout] = INITIALIZATION(p);

    % Run the laser simulation
    Data = DATA_SEQUENCE(Pin, Pulse_duration, p);
    [Ng, Nabs, Nph, Pout] = Laser(p, Data, 0, Vabs);

    % Plot the phase portrait
    figure;
    plot3(Ng, Nabs, Pout, 'LineWidth', 1.5);
    xlabel('Ng');
    ylabel('Nabs');
    zlabel('Pout');
    title('Laser Phase Portrait');
    grid on;
end

% ... (rest of your functions remain unchanged)

% Call the function to generate the phase portrait




function [Pout] = Laser(p,Data,I_bias,Vabs)
    [ Ng , Nabs , Nph , Pout ] = INITIALIZATION( p );

    for cc_t = 1 : p.stab - 1
        % Laser No1     
        Ag = - 1 / p.tg;
        Bg = - p.vg * p.Gamma_g * GAIN( Ng(cc_t) / p.Vg ) * Nph( cc_t ) * p.Vg / ( p.Vg + p.Va ) + I_bias / p.q; 
        Aa = - 1 / p.ta;
        Ba = - p.vg * p.Gamma_a * GAIN( Nabs( cc_t ) / p.Va ) * Nph( cc_t ) * p.Va / ( p.Vg + p.Va ); 
        As = p.vg * p.Gamma_g * GAIN( Ng( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
        Bs = p.bsp* Ng( cc_t ) /p.tg; 

        Ng( cc_t + 1 ) = Ng( cc_t ) * exp( Ag * p.dt ) + Bg / Ag * ( exp( Ag * p.dt ) - 1 );
        Nabs( cc_t + 1 ) = Nabs( cc_t ) * exp( Aa * p.dt ) + Ba / Aa * ( exp( Aa * p.dt ) - 1 );
        Nph( cc_t + 1 ) = Nph( cc_t ) * exp( As * p.dt ) + Bs / As * ( exp( As * p.dt ) - 1 );
        Pout( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph( cc_t + 1 ) / ( p.tph * p.lambda );
    end
    for cc_t = p.stab : p.tot_time
        % Input power
        Pinput = p.eta_input * p.tph * p.lambda * Data( cc_t ) / ( p.h * p.c ); 
    
        % Laser 1    
        Ag = - 1 / p.tg;
        Bg = - p.vg * p.Gamma_g * GAIN( Ng( cc_t ) / p.Vg ) * ( Nph( cc_t ) - Pinput ) * p.Vg / ( p.Vg + p.Va ) + I_bias / p.q; 
        Aa = - 1 / p.ta;
        Ba = - p.vg * p.Gamma_a * GAIN( Nabs( cc_t ) / p.Va ) * Nph( cc_t ) * p.Va / ( p.Vg + p.Va ); 
        As = p.vg * p.Gamma_g * GAIN( Ng( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
        Bs = p.bsp* Ng( cc_t ) /p.tg; 
                        
        Ng( cc_t + 1 ) = Ng( cc_t ) * exp( Ag * p.dt ) + Bg / Ag * ( exp( Ag * p.dt ) - 1 );
        Nabs( cc_t + 1 ) = Nabs( cc_t ) * exp( Aa * p.dt ) + Ba / Aa * ( exp( Aa * p.dt ) - 1 );
        Nph( cc_t + 1 ) = Nph( cc_t ) * exp( As * p.dt ) + Bs / As * ( exp( As * p.dt ) - 1 );
        Pout( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph( cc_t + 1 ) / ( p.tph * p.lambda );
    end
end
function Data = DATA_SEQUENCE( Pin , Pulse_duration , p )
    Data = zeros( 1 , p.tot_time );
    %Data( p.stab + 1 : p.stab + round( Pulse_duration / p.dt ) ) = Pin * ones( 1 , round( Pulse_duration / p.dt ) );
end

function  [ Ng , Nabs , Nph , Pout ] = INITIALIZATION( p )
    Ng = 1e1 * ones( 1 , p.stab );
    Nabs = 1e1 * ones( 1 , p.stab );
    Nph = 1e-10 * ones( 1 , p.stab );
    Pout = 1e-10 * ones( 1 , p.stab );
end

function g = GAIN( n )
    %% From fitting of the Laser Data
    a = 9.571e4;    b= -5.229e6;
    g = a * log( n ) + b;
end
 

function p = constants(Vabs,L,Rga)
    %% global constants
    p.h = 6.626070e-34;
    p.c =  3e8;
    p.q = 1.6021764*1e-19;
    p.T = 300;
    p.kb = 1.380649e-23;
    
    %% materail constants
    p.ng = 3.43;
    p.a = 1000;
    p.a1 = 0.4;
    p.a2 = 0.6;
    p.vg = p.c / p.ng;
    
    %% mirrors
    p.R1 = 1;
    p.R2 = 0.6;
    
    %% dimanesions
    p.H = 98e-9;
    p.W = 4e-6;
    p.Lg = ( 1 - Rga ) * L;
    p.La = Rga * L;
    p.Vg = p.Lg * p.H * p.W;
    p.Va = p.La * p.H * p.W;
    p.V = (p.Lg + p.La)*p.H*p.W;
    p.lambda = 1550e-9;
    p.Gamma_g = 0.09;
    p.Gamma_a = 0.09;
    p.tg = 1e-9;
    recombination_rate = 1/p.tg;
    G0 = 1/100e-12 - 1/p.tg;
    thermal_escape_rate = G0*exp(p.H*p.q*abs(Vabs)/(2*p.W*p.kb*p.T));
    t_sa = 1/(recombination_rate + thermal_escape_rate);
    p.ta =t_sa;
    p.tph = 1 / ( p.c / p.ng * ( p.a - log( p.R1 * p.R2) / ( 2 * ( p.Lg + p.La ) ) ) * ( 1 + p.a1 + p.a2 * (Vabs) ));
    p.Ng0 = 5.11e+23; 
    p.Nabs0 = 5.11e+23; 
    p.Br = 10e-16;
    p.bsp = 1e-4;
    p.eta_c = 0.4;
    p.eta_input = 0.4;
    p.Lc = 2*L/(p.vg*p.tph);
    p.s = 5e-23/p.V;    

    %% time constants 
    p.dt = 1e-12;
    p.stab = round( 10e-9 / p.dt );
    p.tot_time = round( 30e-9 / p.dt );
end

