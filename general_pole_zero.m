%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conor Igoe - Advanced Signal Processing Assignment 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

blue = [0, 0.4470, 0.7410];
red = [0.8500, 0.3250, 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];

upper_z = 1.3;

pos1 = [0.1 0.55 0.4 0.35]; %surface plot
pos2 = [0.6 0.55 0.3 0.35]; %pole zero in z-space
pos3 = [0.6 0.1 0.3 0.35]; %pole zero in s-space
pos4 = [0.1 0.325 0.4 0.125]; %freq response
pos5 = [0.1 0.1 0.4 0.125]; %kernel

nbit = 13;  % coefficient wordlength including sign bit
maxval = 2^(nbit-1) - 1; %Maximum integer value possible using sign magnitude representation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Original, Non-Quantized Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bpFIR = designfilt('bandpassfir','PassbandFrequency1',0.345,...
                              'StopbandFrequency1',0.255,...
                              'PassbandRipple',2,...
                              'StopbandAttenuation1',62,...
                              'PassbandFrequency2',0.655,...
                              'StopbandFrequency2',0.745,...
                              'StopbandAttenuation2',62,...
                              'DesignMethod','equiripple');   %Design for more aggresive, symmetrical filter spec

og_kernel = bpFIR.Coefficients;
og_zeros = roots(og_kernel);

real_time_plot_store = ones(length(og_zeros),1);
best_plot_store = ones(length(og_zeros),1);

s_real_time_plot_store = ones(length(og_zeros),1);
s_best_plot_store = ones(length(og_zeros),1);

best_kernel = ones(1,length(og_kernel));
best_score = 100000;
best_zeros = ones(length(og_kernel),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
subplot('Position',pos2)
[hz,hp,ht] = zplane(og_zeros,zeros(length(og_zeros),1));
set(findobj(hz, 'Type', 'line'), 'Color', 'w'); %Disable default pole zero plotting
title('Pole Zero Plot in z-domain')
xlabel('Re(z)')
ylabel('Im(z)')
grid on
hold on

for index = [1:length(og_zeros)]
    plot(og_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', blue, 'LineWidth',2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Closest Quantized Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxcoeff = max(abs(og_kernel));

q_scale = maxval/maxcoeff;   % Scaling factor
q_kernel = q_scale.*og_kernel;
q_kernel = round(q_kernel);  %Rounded scaled coefficients  (all are now integers)
closest_kernel = q_kernel/q_scale;

closest_zeros = roots(closest_kernel);
closest_zeros = closest_zeros + ones(length(closest_zeros),1).*0.00001i;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate score of closest quantized filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_test_zeros = closest_zeros;
min_weighted_dist = ones(length(og_zeros),1)*100000; %Initialize to large value
remove_index = 0;
for index = 1:length(og_zeros)
    for index2 = 1:length(temp_test_zeros)
        d1 = abs(og_zeros(index));
        d1_dash = abs(temp_test_zeros(index2)); 
        theta1 = angle(og_zeros(index));
        theta1_dash = angle(temp_test_zeros(index2)); 
        weight = abs(log(1 + abs(d1)));
        sq_dist = abs(log(d1_dash) - log(d1))*abs(theta1 - theta1_dash) - abs(log(d1));
        if sq_dist < min_weighted_dist(index)
            min_weighted_dist(index) = sq_dist;
            remove_index = index2;
        end
    end
    temp_test_zeros(remove_index) = [];
end

closest_score = sum(min_weighted_dist);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot initial results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index = 1:length(closest_zeros)
    subplot('Position',pos2)
    plot(closest_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', red, 'LineWidth',2);
    best_plot_store(index) = plot(best_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', purple, 'LineWidth',2);
    subplot('Position',pos3)
    hold on
    plot(log(abs(closest_zeros(index))),angle(closest_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', red, 'LineWidth',2);
    s_best_plot_store(index) = plot(log(abs(best_zeros(index))),angle(best_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', purple, 'LineWidth',2);    
end

for index = 1:length(closest_zeros)
    subplot('Position',pos2)
    real_time_plot_store(index) = plot(closest_zeros(index), '.', 'MarkerSize', 2, 'MarkerEdgeColor', yellow);
    subplot('Position',pos3)
    s_real_time_plot_store(index) = plot(log(abs(closest_zeros(index))),angle(closest_zeros(index)), '.', 'MarkerSize', 2, 'MarkerEdgeColor', yellow);
end

quantized_kernel = closest_kernel;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 1/maxval;
for simulation_index = 1:100
    
    mask = quantized_kernel./quantized_kernel;
    mask(isnan(mask)) = 0;
    mask(1) = 0;
    
    step_kernel = randi([-1,1],1,length(og_kernel)).*step.*mask;
    test_kernel = quantized_kernel + step_kernel;
    test_zeros = roots(test_kernel);
    test_zeros = test_zeros + ones(length(test_zeros),1).*0.00001i;

    temp_test_zeros = test_zeros;
    min_weighted_dist = ones(length(og_zeros),1)*100000;
    remove_index = 0;

    for index = 1:length(og_zeros)
        for index2 = 1:length(temp_test_zeros)
            d1 = abs(og_zeros(index));
            d1_dash = abs(temp_test_zeros(index2)); 
            theta1 = angle(og_zeros(index));
            theta1_dash = angle(temp_test_zeros(index2)); 
            weight = abs(log(1 + abs(d1)));
            sq_dist = abs(log(d1_dash) - log(d1))*abs(theta1 - theta1_dash) - abs(log(d1));
            if sq_dist < min_weighted_dist(index)
                min_weighted_dist(index) = sq_dist;
                remove_index = index2;
            end
        end
        temp_test_zeros(remove_index) = [];
    end

    score = sum(min_weighted_dist);
    plot_color = yellow;
    if score < best_score
        best_score = score;
        best_kernel = test_kernel;
        best_zeros = test_zeros;
        plot_color = purple;
    end

    for index = 1:length(test_zeros)
        subplot('Position',pos2)
        plot(test_zeros(index), '.', 'MarkerSize', 10, 'MarkerEdgeColor', 'black');
        delete(real_time_plot_store(index));
        real_time_plot_store(index) = plot(test_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', plot_color);
        hold on
        
        subplot('Position',pos3)
        title('Pole Zero Plot in s-domain')
        xlabel('Re(s)')
        ylabel('Im(s)')
        grid on
        plot(log(abs(test_zeros(index))),angle(test_zeros(index)), '.', 'MarkerSize', 10, 'MarkerEdgeColor', 'black');
        delete(s_real_time_plot_store(index));
        s_real_time_plot_store(index) = plot(log(abs(test_zeros(index))),angle(test_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', plot_color);
        hold on
        h=get(gca);
        xlim([-1.5,1.5])
        ylim([-3.5,3.5])
        line(h.XLim,[0 0], 'color', 'black');        
        line([0 0],h.YLim, 'color', 'black');        
    end

    for index = 1:length(best_zeros)
        subplot('Position',pos2)
        delete(best_plot_store(index));
        best_plot_store(index) = plot(best_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', purple);
        plot(og_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', blue, 'LineWidth',2);
        plot(closest_zeros(index), '.', 'MarkerSize', 20, 'MarkerEdgeColor', red, 'LineWidth',2);
        hold on
        subplot('Position',pos3)
        delete(s_best_plot_store(index));
        s_best_plot_store(index) = plot(log(abs(best_zeros(index))),angle(best_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', purple);
        plot(log(abs(og_zeros(index))),angle(og_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', blue, 'LineWidth',2);
        plot(log(abs(closest_zeros(index))),angle(closest_zeros(index)), '.', 'MarkerSize', 20, 'MarkerEdgeColor', red, 'LineWidth',2);
        hold on
    end  

    subplot('Position',pos4)
    [w, freq_res_og, freq_res_test] = get_freq_res(og_zeros, test_zeros);
    [w, freq_res_og, freq_res_best] = get_freq_res(og_zeros, best_zeros);
    [w, freq_res_og, freq_res_closest] = get_freq_res(og_zeros, closest_zeros);
    cla
    w = w./pi;
    plot(w,20*log10(abs(freq_res_og))+20*log10(abs(og_kernel(1))),'LineWidth',2);
    hold on
    plot(w,20*log10(abs(freq_res_closest))+20*log10(abs(closest_kernel(1))),'LineWidth',2);
    plot(w,20*log10(abs(freq_res_test))+20*log10(abs(test_kernel(1))),'LineWidth',2);
    plot(w,20*log10(abs(freq_res_best))+20*log10(abs(best_kernel(1))),'LineWidth',2);

    legend('Original','Closest',sprintf('Test (%d)',simulation_index),'Best')
    xlim([0,1])
    ylim([-75,5])
    grid on
    title('Frequency Response')
    xlabel('Normalized Frequency \pi \times rad / sec')
    ylabel('20log(|H(j\omega)|)')

    subplot('Position',pos5)
    cla
    stem(1:length(og_kernel),og_kernel, '.', 'MarkerSize', 12, 'LineWidth', 2);
    title('FIR Filter Kernel')   
    xlabel('Kernel Sample Index')
    ylabel('Kernel Sample Value')
    hold on
    grid on
    stem((1:length(closest_kernel))+0.2,closest_kernel, '.', 'MarkerSize', 12, 'LineWidth', 2);
    stem((1:length(test_kernel))+0.4,test_kernel, '.', 'MarkerSize', 12,'LineWidth', 2);
    stem((1:length(best_kernel))+0.6,best_kernel, '.', 'MarkerSize', 12,'LineWidth', 2);
    legend('Original','Closest',sprintf('Test (%d)',simulation_index),'Best', 'Location', 'southeast')    
    
    subplot('Position',pos1)
    cla
    [x,y] = meshgrid(-1.5 : 0.005: 1.5, -1.5 : 0.005: 1.5);
    s = x + i*y;
    z = ones(length(s));
    for index = 1:length(test_zeros)
        z = z.*(s - test_zeros(index))./s;
    end
    z=abs(z);
    z = z.*abs(test_kernel(1));
    z(z > upper_z) = upper_z + .05;
    surf(x, y, z)   
    title('Pole Zero Topology (z-domain)')
    hold on
    xlim([-1.5,1.5])
    ylim([-1.5,1.5])
    zlim([0,upper_z])   
    
    r = 1;
    theta = 0:0.001:2*pi;
    x = r*cos(theta);
    y = r*sin(theta);
    s = x + i*y;
    z = ones(1,length(s));
    for index = 1:length(test_zeros)
        z = z.*(s - test_zeros(index))./s;
    end
    z=abs(z);
    z = z.*abs(test_kernel(1));
    plot3(x,y,z, 'LineWidth',4, 'Color', 'red');
    for index = 1:length(test_zeros)
        plot3(real(test_zeros(index)),imag(test_zeros(index)), 0.025,'.w', 'MarkerSize', 12)
    end    
    plot(x,y, '--r')
    plot(0,0, '*')
    xlabel('Re(z)')
    ylabel('Im(z)')
    zlabel('|H(z)|')
    axis vis3d
    view(72,45)
    shading interp
    drawnow    
end