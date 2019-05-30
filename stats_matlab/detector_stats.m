% For Octave uncomment the following lines
% pkg load signal

area_list = [];
low_freq_list = [];
high_freq_list = [];
duration_list = [];

% Event extracted from the marked (.mrk) files, those times correspond to fast
% ripples in the signal
event_list_CR_15dB_1 = [4.257324 7.109863 10.126953 12.882324 15.416016 18.293457 20.631348 23.371094 33.691406 36.026367 38.924805 41.888672 47.237793 63.417969 71.552246 74.271484 79.210938 81.780273 90.965820 96.573242 99.051758 105.001465 107.542480 110.891602];
event_list_CR_10dB_2 = [1.00 4.126 8.987 12.518 15.445 18.695 22.164551 27.433105 29.932617 35.999 39.245117 41.329102 43.876465 49.174316 52.139648 54.865723 73.749512 79.445313 81.967773 95.738281 104.247070 107.252969 115.780762]; 
event_list_CR_10dB_3 = [3.4131 6.876465 14.625977 20.131348 26.065918 32.775 36.106445 38.892578 42.418457 47.615234 65.390625 72.854492 75.655762 82.522461 85.334 88.101563 90.402344 93.464355 96.665527 99.223633 102.135742 104.865723 107.899414 119.000000]; 
event_list_CR_10dB_4 = [7.172852 10.022461 12.401855 21.015137 23.554199 26.556641 29.043457 32.575684 36.358887 42.98 45.671387 51.123535 54.360352 57.558594 59.862793 67.495117 70.247559 78.650879 83.891602 95.525879 97.634277 103.762207 106.625000 116.340332]; 
event_list_CR_10dB_5 = [1.000000 3.963867 7.450684 10.694336 14.484863 22.042969 29.776367 33.209473 35.962402 47.086426 49.660645 52.327637 55.044434 61.482422 64.537109 73.010742 81.993652 85.161621 88.656250 102.497070 105.273926 107.748047 113.509277 116.443359]; 
event_list_CR_10dB_6 = [6.814453 17.585938 20.162598 26.109863 32.104980 34.892578 37.496582 40.676758 47.603027 62.627441 64.799316 67.22 70.229004 73.166016 79.695313 85.398438 90.396484 95.905273 99.037109 105.62 108.422852 111.116211 113.956055 116.269531]; 
event_list_CR_10dB_7 = [1.000000 12.813477 15.782715 17.795410 21.356445 24.222168 29.074707 31.748535 41.875488 44.715820 47.964355 51.675293 54.433105 58.115234 60.880859 65.352051 73.489746 80.270996 84.911133 87.471191 98.973633 108.511719 111.939941 116.128906]; 
event_list_CR_10dB_8 = [1.000000 7.343750 13.893066 16.474121 22.165039 24.868 27.583008 30.359375 39.430176 42.352051 45.430176 52.232422 58.599121 66.566895 69.490723 81.396973 84.683105 90.393555 96.772949 99.499023 109.54 111.893555 115.974609 119.000000]; 
event_list_CR_10dB_9 = [15.409668 18.183594 20.500488 26.229492 29.420410 32.928223 36.226563 44.312012 47.258301 52.216309 57.382 60.729 66.856934 73.538086 83.245605 88.821777 91.115234 94.322266 96.477051 102.810547 105.771484 113.605957 116.155273 119.000000]; 
event_list_CR_10dB_10 = [1.000000 3.858398 6.149414 8.846680 11.683105 17.793945 19.450684 24.972168 33.753418 37.008789 42.544922 45.952637 58.564941 61.017090 66.205566 69.566895 72.222656 90.743652 93.189453 97.272461 104.409180 107.476074 110.248047 116.437988];

list_event_list = {event_list_CR_15dB_1 event_list_CR_10dB_2 event_list_CR_10dB_3 event_list_CR_10dB_4 event_list_CR_10dB_5 event_list_CR_10dB_6 event_list_CR_10dB_7 event_list_CR_10dB_8 event_list_CR_10dB_9};

% Get some useful infos for what follows
load ../data/CR-SNR-10dB-2.mat;

fs = sr; % sampling frequency
small_duration = 120; % Duration of the data segment that will be analyzed
t=0:1/fs:small_duration-1/sr; % Time axis

% Now each data file

load ../data/CR5-15dB-1.mat
sig_CR_15dB_1 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-2.mat;
sig_CR_10dB_2 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-3.mat;
sig_CR_10dB_3 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-4.mat;
sig_CR_10dB_4 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-5.mat;
sig_CR_10dB_5 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-6.mat;
sig_CR_10dB_6 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-7.mat;
sig_CR_10dB_7 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-8.mat;
sig_CR_10dB_8 = (data(1:int32(small_duration*sr))/100).'; % actual data
load ../data/CR-SNR-10dB-9.mat;
sig_CR_10dB_9 = (data(1:int32(small_duration*sr))/100).'; % actual data

list_sig = [sig_CR_15dB_1 sig_CR_10dB_2 sig_CR_10dB_3 sig_CR_10dB_4 sig_CR_10dB_5 sig_CR_10dB_6 sig_CR_10dB_7 sig_CR_10dB_8 sig_CR_10dB_9];

for i = 1:8
    s = list_sig(:,i);
    event_list = list_event_list{i};
    fprintf('Signal number %d\n', i+1);

    % Analyse each marked event
    for time = event_list(1:end)
        fprintf('\tEvent at time %d s : ', time);
        
        time_scale = 1.5; % length of the analyses window centered on the EOI
        freq_scale = 1; % controls the frequency range shown in the TF map

        % Compute the time frame around the event
        N0 = int32((time-time_scale/2) * fs);
        N1 = int32((time+time_scale/2) * fs);
        t0 = time - time_scale/2;
        t1 = time + time_scale/2;

        % Select the actual event data
        S = s(N0:N1);
        Time = t(N0:N1);

        % We need a gaussian window to conform to the theory
        size_of_window = floor(fs / 8) + 1 - rem(floor(fs/8),2);
        gaussian_window = tftb_window(size_of_window, 'Gauss');
        [tfr, eater, fr] = tfrstft(S,double(1:N1-N0+1),fs,gaussian_window);

        % Find local minima, using find_zeros.oct(_mex)
        [zero_y zero_x] = find_zeros_mex(abs(tfr));
        zero_x = (zero_x-1) / fs + t0;
        zero_y = (zero_y-1) / fs;
        zero_x = zero_x(find(zero_y <= 1/2));
        zero_y = zero_y(find(zero_y <= 1/2));

        % Use the built-in delaunay triangulation
        T = delaunay (double(zero_x), double(zero_y));

        % Store the triangles in a fashion that makes their printing easy
        X = [zero_x(T(:,1));zero_x(T(:,2));zero_x(T(:,3));zero_x(T(:,1))];
        Y = [ zero_y(T(:,1)); zero_y(T(:,2)); zero_y(T(:,3)); zero_y(T(:,1))];

        % Use a frame to limit the side effect (triangles are way bigger on when
        % they are on the side)
        proportion = 1/10;
        t_inf = t0 + (t1-t0)*proportion;
        t_sup = t1 - (t1-t0)*proportion;
        f_inf = 0.02;
        f_sup = 0.48;

        % We only keep triangles between t_inf and t_sup, and f_inf and f_sup

        Y_clr = Y(:,find(min(X) > t_inf & max(X) < t_sup));
        X_clr = X(:,find(min(X) > t_inf & max(X) < t_sup));
        X_clr = X_clr(:,find(min(Y_clr) > f_inf & max(Y_clr) < f_sup));
        Y_clr = Y_clr(:,find(min(Y_clr) > f_inf & max(Y_clr) < f_sup));

        % Compute the max length of a side in each triangle
        side1 = sqrt((X_clr(1,:)-X_clr(2,:)).^2+(Y_clr(1,:)-Y_clr(2,:)).^2);
        side2 = sqrt((X_clr(2,:)-X_clr(3,:)).^2+(Y_clr(2,:)-Y_clr(3,:)).^2);
        side3 = sqrt((X_clr(3,:)-X_clr(1,:)).^2+(Y_clr(3,:)-Y_clr(1,:)).^2);

        dist_t = max(max(side1, side2), side3);
        
        % Here are our thresholds
        dist_threshold = quantile(dist_t, [0.99]);

        % Plot the time freq representation, with the line showing the 
        % detected EOI
        surface = 10*log10(abs(tfr(1:freq_scale*fs/2,:)));
        imagesc(Time, freq_scale * fs * fr(1:fs/2), surface); axis('xy');
        line([time time], [0 freq_scale*1/2*fs], 'color', [1 1 1] );

        % Draw the frame
        hold on;
        x_frame = [t_inf t_sup t_sup t_inf t_inf];
        y_frame = fs * [f_inf f_inf f_sup f_sup f_inf];
        plot(x_frame, y_frame);
        hold off;

        % Keep only the triangles with one of the longest edges
        X_dist = X_clr(:,find(dist_t>dist_threshold));
        Y_dist = Y_clr(:,find(dist_t>dist_threshold));

        % Compute the connected components of adjacency of those triangles
        [components caract]=find_components_mex(X_dist(1:3,:),fs*Y_dist(1:3,:));

        % Find the components of interest (COI) : components that cross the time
        % of the EOI
        COI = find(caract(1,:) < time & caract(2,:) > time);

        % Find the triangles of interest (TOI) : triangles in some COI
        TOI = find(ismember(components, COI));

        % Draw those triangles
        hold on;
        plot (X_dist(:,TOI), fs*Y_dist(:,TOI), 'b', 'color', [0 0 0]);
        hold off;

        % Select components lying between 150Hz and 620Hz, and reaching at
        % least 250Hz

        HFO = sort(COI);
        HFO = HFO(find(caract(3,HFO)>150&caract(4,HFO)<620&caract(4,HFO)>250));
        
        if (length(HFO) > 0)
            % Get the triangles corresponding to the studed components, assumed 
            % to be a HFO thanks to the marking
            TOI = find(ismember(components, HFO));

            % Compute the area of each triangle
            % TODO : break this one liner into a more elegant computation...

            area_t = abs(0.5*(((X_dist(2,TOI)-X_dist(1,TOI)).*(Y_dist(3,TOI)-Y_dist(1,TOI)))-((X_dist(3,TOI)-X_dist(1,TOI)).*(Y_dist(2,TOI)-Y_dist(1,TOI)))));

            % Compute the total area of the component
            area = sum(area_t);

            % Make some stats about all the HFO components, respectively the
            % area of components, their highest and lowest frequencies and their 
            % duration.

            area_list(end+1) = area;
            high_freq_list(end+1) = max(caract(4,HFO));
            low_freq_list(end+1) = min(caract(3,HFO));
            duration_list(end+1) = max(caract(2,HFO) - caract(1,HFO));

            fprintf('HFO. %d component(s), area : %f \n', length(HFO), area);
        else
            fprintf('nothing interesting... \n');
        end
    end;
end;
