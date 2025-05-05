
function dwdt = ode_function_pattern_generation_barray(~, w, b_array, eta, xr, s, Avw, A, ~)
  
    n_neu = size(Avw, 1);
    dwdt = zeros(n_neu*3, 1);

    for i_neu = 1:n_neu

        b = b_array(i_neu);
        
        Csi = zeros(3, n_neu);
        Csi(1, :) = 2*w(1 + ((1:n_neu)-1)*3)-1;

        teta = Avw(i_neu, :)';

        x = w(1 + (i_neu-1)*3);
        y = w(2 + (i_neu-1)*3);
        z = w(3 + (i_neu-1)*3);

        h = [ -x^3 + b*x^2;
                1-5*x^2; 
               -eta.*s.*xr ];

        rumore_Q = 0;
        dwdt((1:3) + (i_neu-1)*3) = A*[x; y; z] + h + Csi*teta + rumore_Q';

    end
    
end
