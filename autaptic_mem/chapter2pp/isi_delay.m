 
%scanning parameters
Pin = 1e-3;
ISI = 1e-9;            %input(blue line): pulse period...how much time it needs to complete one oscillation
nr_cycles = 5;          %number of cycles, how many square pulses we putting in
dc = 1;               %duty cycle: how much time inside the period the pulse stays on...50% means have period 'high' half 'low'
delay = 60e-9;          %transmision distance between the two Lasers
pulse_delay =10e-9;


%THATS THE MAIN CODE --------------------------------------------------------------------------------------
figure('Position', [0 0 2000 1000])
colors = [[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880];  [0.6350 0.0780 0.1840]; [0.6350 0.0780 0.1840]];

%no point trying 0 and 1 because they dont work at all for memory (spikes not enough power to excite next laser)
for Vabs = 2:5
    Vabs
    I_bias_map1 = [16.19,19.87,24.39,29.94,36.75,45.10,55.37].*1e-3;
    I_bias_map2 = I_bias_map1 + .25e-3;

    I_bias1 = I_bias_map1(Vabs + 1);       
    I_bias2 = I_bias_map2(Vabs + 1);  

    L = 200e-6;                                 
    Rga = 0.1;                                   
    inj = 0.6;                                  
    p = constants(Vabs,L,Rga);
    p.tot_cycles = 700;  

    cycles = 1:10;
    ps = 1e-9:.5e-9:10e-9;

    min_pulse_delay = zeros(length(cycles),1);
    for nr_cycles = cycles
        
        for pulse_delay = ps
            delay = (ISI+pulse_delay)*nr_cycles;    %rough estimation of the delay so that we can count the spikes

            Data = DATA_SEQUENCE(Pin,ISI,nr_cycles,dc,p,pulse_delay);
            [Pout1,Pout2] = lasers(p,Data,I_bias1,I_bias2,delay,nr_cycles,ISI);

            %finding the area where the initial spikes will happen with our input 
            start = p.stab;
            finish = p.stab + round((ISI+pulse_delay)/p.dt + 1)*nr_cycles;

            [peaks1_init] = findpeaks(Pout1(start:finish),'MinPeakProminence',0.025,'MinPeakHeight',0.02);
            [peaks2,loc2,~,~] = findpeaks(Pout2,'MinPeakProminence',0.025,'MinPeakHeight',0.02);
            
            %if the number of spikes is not a integer multiplier of our initial number of spikes there is no way to procceed
            %here 3 is choosen as integer but any other number will work as long as we are still in p.total_cycles
            if(length(loc2) < 3*nr_cycles + 1)
                continue;
            end

            %between the last spike of Slave laser and the first of the next batch we should have all the spikes from ML
            start = loc2(3*nr_cycles);
            finish = loc2(3*nr_cycles + 1);
            [peaks1_mid] = findpeaks(Pout1(start:finish),'MinPeakProminence',0.025,'MinPeakHeight',0.02);

            %if the initial spikes are the same after some iterrations (3 choosen here) then the memory works fine and dont need more distance between pulses
            if(length(peaks1_init) == length(peaks1_mid))
                min_pulse_delay(nr_cycles) = pulse_delay;
                break;
            end
        end
    end
    hold on;
    plot(cycles,min_pulse_delay*1e9,'LineWidth' , 3 , 'Color',colors(Vabs+1,:));

    path = ['writting/ps_dist/pulse_distance.fig'];
    savefig(path)
    
end

xlabel('spike count','FontSize' , 20);
ylabel('minimum distance (ps)','FontSize' , 20);
title('minimum pulse distance for pulse regeneration');
legend('Vabs = 0 ', 'Vabs = 1','Vabs = 2', 'Vabs = 3','Vabs = 4','Vabs = 5')

%figure

%% Plot Data
%t = p.dt : p.dt : p.stab * p.dt + ISI * ( nr_cycles + p.tot_cycles );
%plot( Data * 1e3 , 'LineWidth' , 3 )
%hold on
%plot( Pout1 * 1e3 , 'LineWidth' , 3 )
%plot( Pout2 * 1e3 , 'LineWidth' , 3 ,'Color', 'black')
%xlabel( 'time (ns)' , 'FontSize' , 20 )
%ylabel( 'Output Power (mW)' , 'FontSize' , 20)
%legend( 'Input' , 'Neuron 1' , 'Neuron 2' , 'FontSize' , 20 )
%title( [ 'Pin=' num2str( Pin * 1e3 ) 'DC=' num2str( dc ) 'Period=' num2str( ISI * 1e9 ) ] )




function [Pout1,Pout2] = lasers(p,Data,I_bias1,I_bias2,delay,Nr_cycles,ISI)

    [ Ng1 , Nabs1 , Nph1 , Pout1 ] = INITIALIZATION( p );
    [ Ng2 , Nabs2 , Nph2 , Pout2 ] = INITIALIZATION( p );
    
    
    
    %Stabilization for both lasers
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

    %Lasers
    p.Tdelay = round( delay / p.dt );
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
    
    
end
function Data = DATA_SEQUENCE( Pin , ISI , Nr_cycles , DC , p,pulse_delay )
    Data = zeros( 1 , p.stab + round( ISI / p.dt ) * ( Nr_cycles + p.tot_cycles ) );

    for cc_cycle = 1 : Nr_cycles
        start = p.stab + (cc_cycle -1)*round(pulse_delay/p.dt) +( cc_cycle - 1 ) * round( ISI / p.dt ) + 1 ;
        finish = start + round( DC * ISI/ p.dt );

        Data( start : finish ) = Pin * ones();
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
    p.stab = round( 300e-9 / p.dt );
end


