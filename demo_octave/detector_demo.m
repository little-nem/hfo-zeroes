pkg load signal

load ../data/CR-SNR-10dB-2.mat;

% Load the data
fs = sr; % sampling frequency
small_duration = 100; % Duration of the data segment that will be analyzed

t = 0 : 1/fs : small_duration-1/sr; % Time axis
s=data(1:int32(small_duration*sr)); % Actual data

% Add some artefacts, for demonstration purpose, at t=14s and t=18s
spike1 = 3*normpdf(t, 14, 0.004);
spike2 = 0.5*normpdf(t, 18, 0.001);
s = s + spike1 + spike2;

% Transpose s for the rest of the execution
s=s.'; 

% Definitions of the bands we are interested in
ripple_band = [80 250];
fast_ripple_band = [250 500];

% Select the event of interests in some given band, here the fast ripple one
event_list = event_of_interest(s, fs, fast_ripple_band, 4, 0.250);

% Check each event (we reject the first one to avoid potential edge effects)
for time = event_list(2:end)
    fprintf('Event at time %d s : ', time);

    time_scale = 1.5; % length of the analysis window centered on the EOI
    freq_scale = 1; % controls the frequency range shown in the TF map

    % Compute the time frame
    N0 = int32((time-time_scale/2) * fs);
    N1 = int32((time+time_scale/2) * fs);
    t0 = time - time_scale/2;
    t1 = time + time_scale/2;

    % Select the event data
    S = s(N0:N1);
    Time = t(N0:N1);

    % We need a gaussian window to conform to the theory
    size_of_window = floor(fs / 8) + 1 - rem(floor(fs/8),2);
    gaussian_window = tftb_window(size_of_window, 'Gauss');
    [tfr, eater, fr] = tfrstft(S,double(1:N1-N0+1),fs,gaussian_window);
  
    % Find local minima, using find_zeros.oct
    [zero_y zero_x] = find_zeros_oct(abs(tfr));

    zero_x = (zero_x-1) / fs + t0;
    zero_y = (zero_y-1) / fs;
    zero_x = zero_x(find(zero_y <= 1/2));
    zero_y = zero_y(find(zero_y <= 1/2));

    % Use the built-in delaunay triangulation
    T = delaunay (double(zero_x), double(zero_y));

    % Store the triangles in a fashion that makes their printing easy
    X = [zero_x(T(:,1));zero_x(T(:,2));zero_x(T(:,3));zero_x(T(:,1))];
    Y = [ zero_y(T(:,1)); zero_y(T(:,2)); zero_y(T(:,3)); zero_y(T(:,1))];

    % Use a frame to consider only triangle in the center of the screen
    proportion = 1/10;
    t_inf = t0 + (t1-t0)*proportion;
    t_sup = t1 - (t1-t0)*proportion;
    f_inf = 0.02;
    f_sup = 0.48;

    % We only keep triangles between t_inf and t_sup and between f_inf and f_sup

    Y_clear = Y(:,find(min(X) > t_inf & max(X) < t_sup));
    X_clear = X(:,find(min(X) > t_inf & max(X) < t_sup));
    X_clear = X_clear(:,find(min(Y_clear) > f_inf & max(Y_clear) < f_sup));
    Y_clear = Y_clear(:,find(min(Y_clear) > f_inf & max(Y_clear) < f_sup));

    % Compute the length of each side, and keep the longest for each triangle
    side1 = sqrt((X_clear(1,:)-X_clear(2,:)).^2+(Y_clear(1,:)-Y_clear(2,:)).^2);
    side2 = sqrt((X_clear(2,:)-X_clear(3,:)).^2+(Y_clear(2,:)-Y_clear(3,:)).^2);
    side3 = sqrt((X_clear(3,:)-X_clear(1,:)).^2+(Y_clear(3,:)-Y_clear(1,:)).^2);
    dist_t = max(max(side1,side2),side3);
    
    % Here is our threshold
    dist_threshold = quantile(dist_t, [0.99]);

    % Plot the time freq map
    surface = 10*log10(abs(tfr(1:freq_scale*fs/2,:)));
    imagesc(Time,freq_scale*fs*fr(1:fs/2), surface); axis('xy');

    % Draw a line showing the EOI time
    line([time time], [0 freq_scale*1/2*fs], 'color', [1 1 1] );

    % Draw the frame
    hold on;
    plot([t_inf t_sup t_sup t_inf t_inf], fs*[f_inf f_inf f_sup f_sup f_inf]);
    hold off;

    % Keep only the triangles with one of the longest edges
    X_dist = X_clear(:,find(dist_t>dist_threshold));
    Y_dist = Y_clear(:,find(dist_t>dist_threshold));
    
    % Use find_components_oct to compute components and their caracteristics
    % Remember that find_components_oct only needs the first three lines of our
    % representation of the triangles (the extra lines is useful in Matlab to 
    % plot the triangles with a single command, but does not carry any extra 
    % information).

    [components caract] = find_components_oct(X_dist(1:3,:),fs * Y_dist(1:3,:));

    % Find the components of interest (COI) : components that cross the time
    % of the EOI, caract(1,i) being the beginning time of the component i, and 
    % caract(2,i) being its ending time.

    COI = find(caract(1,:) < time & caract(2,:) > time);

    % Keep triangles in some COI to stay focused on the EOI, and not the whole 
    % time frame
    TOI = find(ismember(components, COI));

    hold on;
    plot (X_dist(:,TOI), fs*Y_dist(:,TOI), 'b', 'color', [0 0 0]);
    hold off;

    % Example of criterion we can think about once the caracteristics of the
    % components have been extracted : we look for a component that lies between
    % 200Hz and 650Hz. It's just a proof of feasibility here, and we might want
    % to refine such criterions.

    if (length(find(caract(3,COI) > 200 & caract(4,COI) < 650)) > 0)
       fprintf('HFO: fast ripple \n');
    else
        fprintf('nothing interesting... \n');
    end;

    pause;
end;
