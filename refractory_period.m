figure('Position', [0 0 2560 1600]);    

for Vabs = 0:2:6
    %params
    L = 200e-6;
    Rga = 0.1;
    inj = 0.6;
    I_bias_map = [16.19,19.87,24.39,29.94,36.75,45.10,55.37].*1e-3;  %I_bias for operating as neuron for each Vabss from 0-5
    I_bias = I_bias_map(Vabs+1);
    p = constants(Vabs,L,Rga);


    Pin = 2e-3;
    Pulse_duration = 1e-9;



    ds = 1e-9:5e-10:9e-9;
    pks = zeros(length(ds),3);
    abs_refract = 0;
    rel_refract = 0;

    counter = 1;
    for distance = ds 
        
        Data = DATA_SEQUENCE_2PULSE( Pin , Pulse_duration , p, distance );
        Pout = Laser(p,Data,I_bias,Vabs);

        [peaks,~,~,prom] = findpeaks(Pout(10000:end),'MinPeakProminence',0.025,'MinPeakHeight',0.02);
        pks(counter,1:length(peaks)) = peaks;
        counter = counter + 1;
        if(abs_refract == 0 && length(peaks) == 2)
            abs_refract = distance;
        end
        if(rel_refract == 0 && length(peaks) == 2)
            if(abs(peaks(1) - peaks(2) < 0.000001))
                rel_refract = distance;
            end
        end

    end


    subplot(4,1,Vabs/2 + 1);
    legend( ['absolute refractory period = ' num2str(abs_refract/p.dt) 'S' 'relative refractory period = ' num2str(rel_refract/p.dt) 'S'],'Location', 'northwest' , 'FontSize' , 20 )
    abs_refract
    rel_refract
    for i = 1:2
        plot(ds/p.dt*1e-3,smooth(pks(:,i).')*1e3,'LineWidth' , 3);
        hold on;
    end
    title(['Vabs = ' num2str(Vabs)])
    ylabel( 'Peak Power (mW)' , 'FontSize' , 20)
    xlabel( 'Distance (ns)' , 'FontSize' , 20 )

    
    
end
exportgraphics(gcf,['writting/chapter1/refractory/ref.png'])

%% Plot Data
t = p.dt : p.dt : p.tot_time * p.dt;
figure
plot( Data * 1e3 , 'LineWidth' , 3 )
hold on
plot( Pout * 1e3 , 'LineWidth' , 3 )
xlabel( 'time (ns)' , 'FontSize' , 20 )
ylabel( 'Output Power (mW)' , 'FontSize' , 20)
legend( 'Input' , 'Neuron' , 'FontSize' , 20 )
title( [ 'Pin=' num2str( Pin * 1e3 ) 'Pulse Duration=' num2str( Pulse_duration * 1e9 ) ] )



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
    Data( p.stab + 1 : p.stab + round( Pulse_duration / p.dt ) ) = Pin * ones( 1 , round( Pulse_duration / p.dt ) );
end
function Data = DATA_SEQUENCE_2PULSE (Pin, Pulse_duration,p,distance)
    Data = zeros( 1 , p.tot_time );
    start1 = p.stab + 1 ;
    end1 = start1- 1 + round( Pulse_duration / p.dt );
    %create the first pulse
    Data(start1 : end1) = Pin * ones( 1 , round( Pulse_duration / p.dt ) );
    start2 = end1 + round(distance/p.dt) + 1; 
    end2 = start2 - 1 + round(Pulse_duration/p.dt);
    %create the second pulse
    Data (start2 : end2) =  Pin * ones( 1 , round( Pulse_duration / p.dt ) );
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
    p.eta_input = 0.4;
    p.Lc = 2*L/(p.vg*p.tph);
    p.s = 5e-23/p.V;    

    %% time constants 
    p.dt = 1e-12;
    p.stab = round( 10e-9 / p.dt );
    p.tot_time = round( 30e-9 / p.dt );
end