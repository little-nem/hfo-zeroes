function event_list = event_of_interest(s, sampling_fr, band, k, time_threshold)
% event_list = event_of_interest(s, fs, band, k, time_threshold);
%   Returns a list of event of interest (the times at which each event happens)
%   Method : filter the signal, compute its Hilbert enveloppe and grab events
%   above a threshold given by k times the standard deviation of the enveloppe
%   amplitude, the event list is then made sparse by removing events close to
%   eachother in time, thus corresponding to a single physiological event
%
% inputs
%   s               = Input Signal
%   sampling_fr     = Sampling frequency of the input signal
%   band            = Array of the two frequencies that delimit the band we
%                       are interested in. (eg. [250 500] for fast ripples)
%   k               = Factor in front of the threshold for envelop energy
%                       that is k * std((filtered(s)))
%   time_threshold  = Minimum time between two distincts events (events too
%                       close are merged
%
% output
%   event_list      = Array of time at which event of interest occurs 
%
% Calls:
%   filtfilt
%   hilbert

    if nargin < 3
        error("event_list needs at least 3 arguments...\n");
    elseif nargin == 3
        k = 4;
        time_threshold = 0.250; 
    elseif nargin == 4
        time_threshold = 0.250;
    endif
    
    band_down = band(1);
    band_up = band(2);

    half_fs = sampling_fr/2;

    filtered_s = filtfilt(fir1(64, [band_down/half_fs band_up/half_fs]), 1, s);

    % compute standard deviation of the filtered signal
    sigma = std(filtered_s);

    % find values superior to the threshold
    threshold = k*sigma;
    idx = find(abs(hilbert(filtered_s)) > threshold);

    % Clean the idx : make it sparse (keep only the beginning of the event, the 
    % maximum duration of a single event being time_threshold)
    event_list = [idx(1)/sampling_fr];
    for x = (idx/sampling_fr).'
        if(x-event_list(end) > time_threshold)
            event_list(end+1) = x;
        endif
    endfor
end;
