
clear
close all

start_folder = pwd;
rng(5)

%% time vector 
tspan = [0, 500];
dt = 0.01;
n_sample = tspan(2)/dt;
tspan = linspace(tspan(1), tspan(2), n_sample);

%% neurons
n_neuron_array = 2;
grado_di_em_nonem_array = 0.25;

%% successful rate
SR_matrix   = zeros(length(n_neuron_array) + 1, length(grado_di_em_nonem_array) + 1);
SR_matrix(1, 2:end) = grado_di_em_nonem_array;
SR_matrix(2:end, 1) = n_neuron_array;

%% iteration parameters
n_iterazione = 0;
n_iterazioni_tot = length(n_neuron_array).*length(grado_di_em_nonem_array);
i_neu = 1;

%% percentages of emulation/non-emulation; bursting/spiking
perc_em     = 0;
perc_non_em = 1-perc_em;
perc_burst  = 1;
perc_spike  = 1-perc_burst;

%% iterations
for n_neuron = n_neuron_array

    i_neu = i_neu + 1;
    
    for i_grado_di_em_nonem = 1:length(grado_di_em_nonem_array)
    
        grado_di_em_nonem = grado_di_em_nonem_array(i_grado_di_em_nonem);
        Avw = zeros(n_neuron);
        for i_colonna = 1:(n_neuron-1)
            Avw(randi(i_colonna), i_colonna + 1) = 1;
        end
        n_relazioni = sum(Avw(:) == 1);
        n_rel_em = round(perc_em * n_relazioni);
        n_rel_non_em = n_relazioni - n_rel_em;
        indici_relazioni = Avw == 1;
        vettore_relazioni_da_modificare = Avw(indici_relazioni);
        indici_relazioni_aggiornate = randperm(n_relazioni, n_rel_non_em);
        vettore_relazioni_aggiornate = vettore_relazioni_da_modificare;
        vettore_relazioni_aggiornate(indici_relazioni_aggiornate) = -1;
        Avw(indici_relazioni) = vettore_relazioni_aggiornate;
        
        Avw = Avw + Avw';
        Avw = Avw .* grado_di_em_nonem;

        if size(Avw, 1) ~= n_neuron
            disp("Change matrix Avw")
            return
        end

        %% b array
        n_neu_burst = ceil(n_neuron.*perc_burst);
        n_neu_spike = n_neuron - n_neu_burst;

        b_burst = 2.5 .* ones(1, n_neu_burst);
        b_spike = 3   .* ones(1, n_neu_spike);
        b_array = [b_burst, b_spike];
        
        %% HR parameters
        eta = 0.01;  
        xr = -1;  
        s = 4;
        
        %% HR system
        A = [0,         1,     -1;
             0,         -1,     0;
             eta.*s,    0,     -eta];
        
        n_eq_diff = size(A,1);
        
        %% initial condition of the state
        initial_conditions = rand(n_neuron*n_eq_diff, 1);
        
        %% guess parameters (initial condition of estimated parameters)
        Avw_guess = rand(n_neuron); 
        
        %% output matrix 
        C = [1, 0, 0];
        n_uscite = size(C,1);
                
        %% Q and R values        
        Q = 2.*0.001.*eye(3);
        R = 0.01;
        
        %% observer parameters
        gamma = 1.*eye(n_neuron);
        sigma = 1;
        
        %% pattern generation
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        [t_gen, real_w] = ode45(@(t_gen_i, w) ode_function_pattern_generation_barray(t_gen_i, w, b_array, eta, xr, s, Avw, A, Q), tspan, initial_conditions);
        
        %% plot results
        real_w_observable = real_w(:, 1 + 3*((1:n_neuron)-1));
        
        figure
        title('True values only', 'FontSize', 18)
        hold on  
        for i_neu = 1:n_neuron
            % plot(t, real_w(:, 1 + 3*(i_neu-1)) )
            plot(t_gen, real_w_observable(:, i_neu), 'LineWidth', 1)
        end

        box off
        set(gca, 'TickDir', 'out')
        set(gca, 'LineWidth', 1.2)
        legend('Neuron 1', 'Neuron 2')
        title('Pattern Generation')
        ylabel('Intracellular Membrane Potential')
        xlabel('Time (samples)')
        ylim([-3, 3])

        n_iterazione = n_iterazione + 1;
        disp(strcat("iteration ", string(n_iterazione), " of ", string(n_iterazioni_tot)))


    end

end

cd(start_folder)

