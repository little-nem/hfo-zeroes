function p = normpdf(x,m,s)
% Returns the classic normalpdf, to avoid the nan package which provides
% low performance and shadows some stat functions, at least in Octave...
    if nargin==1,
            m=0;s=1;
    elseif nargin==2,
            s=1;
    end;        

    z = (x-m)./s; % if this line causes an error, input arguments do not fit. 
    SQ2PI = 2.5066282746310005024157652848110;
    p = exp(-z.^2/2)./(s*SQ2PI);
