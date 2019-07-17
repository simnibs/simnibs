function [EC] = spm_ECdensity(STAT,t,df)
% Returns the Euler characteristic (EC) density
% FORMAT function [EC] = spm_ECdensity(STAT,t,df)
%__________________________________________________________________________
%
% Reference : Worsley KJ et al (1996), Hum Brain Mapp. 4:58-73
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ECdensity.m 5544 2013-06-12 11:01:49Z guillaume $


% EC densities
%--------------------------------------------------------------------------
t = t(:)';
if      STAT == 'Z'

    % Gaussian Field
    %----------------------------------------------------------------------
    a       = 4*log(2);
    b       = exp(-t.^2/2);

    EC(1,:) = 1 - spm_Ncdf(t);
    EC(2,:) = a^(1/2)/(2*pi)*b;
    EC(3,:) = a/((2*pi)^(3/2))*b.*t;
    EC(4,:) = a^(3/2)/((2*pi)^2)*b.*(t.^2 - 1);

elseif  STAT == 'T'

    % T - Field
    %----------------------------------------------------------------------
    v       = df(2);
    a       = 4*log(2);
    b       = exp(gammaln((v+1)/2) - gammaln(v/2));
    c       = (1+t.^2/v).^((1-v)/2);

    EC(1,:) = 1 - spm_Tcdf(t,v);
    EC(2,:) = a^(1/2)/(2*pi)*c;
    EC(3,:) = a/((2*pi)^(3/2))*c.*t/((v/2)^(1/2))*b;
    EC(4,:) = a^(3/2)/((2*pi)^2)*c.*((v-1)*(t.^2)/v - 1);

elseif  STAT == 'X'

    % X - Field
    %----------------------------------------------------------------------
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = t.^(1/2*(v - 1)).*exp(-t/2-gammaln(v/2))/2^((v-2)/2);

    EC(1,:) = 1 - spm_Xcdf(t,v);
    EC(2,:) = a^(1/2)*b;
    EC(3,:) = a*b.*(t-(v-1));
    EC(4,:) = a^(3/2)*b.*(t.^2-(2*v-1)*t+(v-1)*(v-2));

elseif  STAT == 'F'

    % F Field
    %----------------------------------------------------------------------
    k       = df(1);
    v       = df(2);
    a       = (4*log(2))/(2*pi);
    b       = gammaln(v/2) + gammaln(k/2);

    EC(1,:) = 1 - spm_Fcdf(t,df);
    EC(2,:) = a^(1/2)*exp(gammaln((v+k-1)/2)-b)*2^(1/2)...
              *(k*t/v).^(1/2*(k-1)).*(1+k*t/v).^(-1/2*(v+k-2));
    EC(3,:) = a*exp(gammaln((v+k-2)/2)-b)*(k*t/v).^(1/2*(k-2))...
              .*(1+k*t/v).^(-1/2*(v+k-2)).*((v-1)*k*t/v-(k-1));
    EC(4,:) = a^(3/2)*exp(gammaln((v+k-3)/2)-b)...
              *2^(-1/2)*(k*t/v).^(1/2*(k-3)).*(1+k*t/v).^(-1/2*(v+k-2))...
              .*((v-1)*(v-2)*(k*t/v).^2-(2*v*k-v-k-1)*(k*t/v)+(k-1)*(k-2));
end
