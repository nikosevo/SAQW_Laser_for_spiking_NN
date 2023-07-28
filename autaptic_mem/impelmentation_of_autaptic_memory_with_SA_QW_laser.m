%% set laser parameters
I_bias1 = 24.45e-3;                          % bias of the gain section 1st laser
I_bias2 = 24.7e-3;                             % bias of the gain section 2nd laser
Vabs = 2;                                    % Absorber's reverse bias 
L = 200e-6;                                  % Length of the laser
Rga = 0.1;                                   % Absorber to laser's length ratio
inj = 0.6;                                   % injection strength
p = constants( Vabs , L , Rga );

%% set scanning parameters
Pin_min = 5e-3;       step_Pin = 0.1e-3;     Pin_max = 5e-3;            % input power of data
ISI_min = 3e-9;       step_ISI = 1e-9;       ISI_max = 3e-9;            % Period of the square pulses
Nr_cycles_min = 3;    step_Nr_cycles = 2;    Nr_cycles_max = 3;         % Number of input square pulses
DC_min = 0.5;         step_DC = 0.1;         DC_max = 0.5;              % Duty Cycle of the pulses
Delay_min = 100e-9;    step_Delay = 1e-9;     Delay_max = 100e-9;         % Transmission distance between the two lasers 
p.tot_cycles = 200;                                                     % Number of simulated pulses
           
%% Initialize laser
[ Ng1 , Nabs1 , Nph1 , Pout1 ] = INITIALIZATION( p );
[ Ng2 , Nabs2 , Nph2 , Pout2 ] = INITIALIZATION( p );

%% Stabilization of the output of the two lasers
for cc_t = 1 : p.stab - 1
    % Laser No1     
    Ag1 = - 1 / p.tg;
    Bg1 = - p.vg * p.Gamma_g * GAIN( Ng1(cc_t) / p.Vg ) * Nph1( cc_t ) * p.Vg / ( p.Vg + p.Va ) + I_bias1 / p.q; 
    Aa1 = - 1 / p.ta;
    Ba1 = - p.vg * p.Gamma_a * GAIN( Nabs1( cc_t ) / p.Va ) * Nph1( cc_t ) * p.Va / ( p.Vg + p.Va ); 
    As1 = p.vg * p.Gamma_g * GAIN( Ng1( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs1( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
    Bs1 = p.bsp* Ng1( cc_t ) /p.tg; 

    Ng1( cc_t + 1 ) = Ng1( cc_t ) * exp( Ag1 * p.dt ) + Bg1 / Ag1 * ( exp( Ag1 * p.dt ) - 1 );
    Nabs1( cc_t + 1 ) = Nabs1( cc_t ) * exp( Aa1 * p.dt ) + Ba1 / Aa1 * ( exp( Aa1 * p.dt ) - 1 );
    Nph1( cc_t + 1 ) = Nph1( cc_t ) * exp( As1 * p.dt ) + Bs1 / As1 * ( exp( As1 * p.dt ) - 1 );
    Pout1( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph1( cc_t + 1 ) / ( p.tph * p.lambda );

    % Laser No2     
    Ag2 = - 1 / p.tg;
    Bg2 = - p.vg * p.Gamma_g * GAIN( Ng2(cc_t) / p.Vg ) * Nph2( cc_t ) * p.Vg / ( p.Vg + p.Va ) + I_bias2 / p.q; 
    Aa2 = - 1 / p.ta;
    Ba2 = - p.vg * p.Gamma_a * GAIN( Nabs2( cc_t ) / p.Va ) * Nph2( cc_t ) * p.Va / ( p.Vg + p.Va ); 
    As2 = p.vg * p.Gamma_g * GAIN( Ng2( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs2( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
    Bs2 = p.bsp* Ng2( cc_t ) /p.tg; 
    
    Ng2( cc_t + 1 ) = Ng2( cc_t ) * exp( Ag2 * p.dt ) + Bg2 / Ag2 * ( exp( Ag2 * p.dt ) - 1 );
    Nabs2( cc_t + 1 ) = Nabs2( cc_t ) * exp( Aa2 * p.dt ) + Ba2 / Aa2 * ( exp( Aa2 * p.dt ) - 1 );
    Nph2( cc_t + 1 ) = Nph2( cc_t ) * exp( As2 * p.dt ) + Bs2 / As2 * ( exp( As2 * p.dt ) - 1 );
    Pout2( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph2( cc_t + 1 ) / ( p.tph * p.lambda );
end

for Pin = Pin_min : step_Pin : Pin_max
    for ISI = ISI_min : step_ISI : ISI_max
        for Nr_cycles = Nr_cycles_min : step_Nr_cycles : Nr_cycles_max
            for DC = DC_min : step_DC : DC_max
                %% Create Input Data stream
                Data = DATA_SEQUENCE( Pin , ISI , Nr_cycles , DC , p );
                   
                for Delay = Delay_min : step_Delay : Delay_max
                    p.Tdelay = round( Delay / p.dt ); 
                    for cc_t = p.stab : p.stab + ( Nr_cycles + p.tot_cycles ) *  round( ISI / p.dt ) - 1
                        % Input power
                        Pin1_to_2 = p.eta_c * p.Gamma_g * p.h * p.c * Nph1( cc_t - p.Tdelay ) / ( p.tph * p.lambda );
                        Pin1_to_2 = p.eta_input * p.tph * p.lambda * Pin1_to_2 / ( p.h * p.c ); 
                        Pin2_to_1 = p.eta_c * p.Gamma_g * p.h * p.c * Nph2( cc_t - p.Tdelay ) / ( p.tph * p.lambda );
                        Pin2_to_1 = p.eta_input * p.tph * p.lambda * Pin2_to_1 / ( p.h * p.c ); 
                        Pinput = p.eta_input * p.tph * p.lambda * Data( cc_t ) / ( p.h * p.c ); 
    
                        % Laser 1    
                        Ag1 = - 1 / p.tg;
                        Bg1 = - p.vg * p.Gamma_g * GAIN( Ng1( cc_t ) / p.Vg ) * ( Nph1( cc_t ) - Pinput - Pin2_to_1 ) * p.Vg / ( p.Vg + p.Va ) + I_bias1 / p.q; 
                        Aa1 = - 1 / p.ta;
                        Ba1 = - p.vg * p.Gamma_a * GAIN( Nabs1( cc_t ) / p.Va ) * Nph1( cc_t ) * p.Va / ( p.Vg + p.Va ); 
                        As1 = p.vg * p.Gamma_g * GAIN( Ng1( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs1( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
                        Bs1 = p.bsp* Ng1( cc_t ) /p.tg; 
                                       
                        Ng1( cc_t + 1 ) = Ng1( cc_t ) * exp( Ag1 * p.dt ) + Bg1 / Ag1 * ( exp( Ag1 * p.dt ) - 1 );
                        Nabs1( cc_t + 1 ) = Nabs1( cc_t ) * exp( Aa1 * p.dt ) + Ba1 / Aa1 * ( exp( Aa1 * p.dt ) - 1 );
                        Nph1( cc_t + 1 ) = Nph1( cc_t ) * exp( As1 * p.dt ) + Bs1 / As1 * ( exp( As1 * p.dt ) - 1 );
                        Pout1( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph1( cc_t + 1 ) / ( p.tph * p.lambda );
    
                        % Laser 2    
                        Ag2 = - 1 / p.tg;
                        Bg2 = - p.vg * p.Gamma_g * GAIN( Ng2( cc_t ) / p.Vg ) * ( Nph2( cc_t ) - Pin1_to_2 ) * p.Vg / ( p.Vg + p.Va ) + I_bias2 / p.q; 
                        Aa2 = - 1 / p.ta;
                        Ba2 = - p.vg * p.Gamma_a * GAIN( Nabs2( cc_t ) / p.Va ) * Nph2( cc_t ) * p.Va / ( p.Vg + p.Va ); 
                        As2 = p.vg * p.Gamma_g * GAIN( Ng2( cc_t ) / p.Vg ) * ( p.Vg / ( p.Vg + p.Va ) ) + p.vg * p.Gamma_a * GAIN( Nabs2( cc_t ) / p.Va ) * ( p.Va / ( p.Vg + p.Va ) ) - 1 / p.tph;
                        Bs2 = p.bsp* Ng2( cc_t ) /p.tg; 
                                       
                        Ng2( cc_t + 1 ) = Ng2( cc_t ) * exp( Ag2 * p.dt ) + Bg2 / Ag2 * ( exp( Ag2 * p.dt ) - 1 );
                        Nabs2( cc_t + 1 ) = Nabs2( cc_t ) * exp( Aa2 * p.dt ) + Ba2 / Aa2 * ( exp( Aa2 * p.dt ) - 1 );
                        Nph2( cc_t + 1 ) = Nph2( cc_t ) * exp( As2 * p.dt ) + Bs2 / As2 * ( exp( As2 * p.dt ) - 1 );  
                        Pout2( cc_t + 1 ) = p.eta_c * p.Gamma_g * p.h * p.c * Nph2( cc_t + 1 ) / ( p.tph * p.lambda );
                    end
    
                    %% Plot Data
                    t = p.dt : p.dt : p.stab * p.dt + ISI * ( Nr_cycles + p.tot_cycles );
                    figure
                    plot( t , Data * 1e3 , 'LineWidth' , 3 )
                    hold on
                    plot( t , Pout1 * 1e3 , 'LineWidth' , 3 )
                    plot( t , Pout2 * 1e3 , 'LineWidth' , 3 )
                    xlabel( 'time (ns)' , 'FontSize' , 20 )
                    ylabel( 'Output Power (mW)' , 'FontSize' , 20)
                    legend( 'Input' , 'Neuron 1' , 'Neuron 2' , 'FontSize' , 20 )
                    title( [ 'Pin=' num2str( Pin * 1e3 ) 'DC=' num2str( DC ) 'Period=' num2str( ISI * 1e9 ) ] )
                end
            end
        end
    end
end

function Data = DATA_SEQUENCE( Pin , ISI , Nr_cycles , DC , p )
    Data = zeros( 1 , p.stab + round( ISI / p.dt ) * ( Nr_cycles + p.tot_cycles ) );

    for cc_cycle = 1 : Nr_cycles
        Data( p.stab + ( cc_cycle - 1 ) * round( ISI / p.dt ) + 1 : p.stab + ( cc_cycle - 1 ) * round( ISI / p.dt ) + round( DC * ISI / p.dt ) ) = Pin * ones( 1 , round( DC * ISI / p.dt ) );
    end
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

function p = constants( Vabs , L , Rga )
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
    p.eta_input = 0.8;
    p.s = 5e-23/p.V;    

    %% time constants 
    p.dt = 1e-12;
    p.stab = round( 20e-9 / p.dt );

end