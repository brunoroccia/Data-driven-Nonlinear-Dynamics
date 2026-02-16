function F = MyForce(t, DATA)
%MyForce External forcing definitions for the Duffing example.
%
%   F = MyForce(t, DATA)
%
% Supports scalar or vector t.
%
% DATA fields used:
%   F    : amplitude
%   W    : angular frequency
%   Flag : selects forcing type
%          0 -> cosine:             F = Amp*cos(W*t)
%          1 -> sine with ramp-on:  F = Amp*sin(W*t)*H(t)
%          2 -> sigmoid ramp:       F = Amp/(1+exp(-k*(t-A))) for t>=A
%          3 -> C^2 smooth step:    transitions from 0 to Amp over [A,C]

Amp = DATA.F;
W   = DATA.W;

switch DATA.Flag
    case 0
        % Pure harmonic forcing
        F = Amp * cos(W*t);

    case 1
        % Sine forcing with C^2 ramp-on over one period
        % Ramp starts at t0 and finishes at t1=t0+2*pi/W
        t0 = 0;
        delta = 2*pi/W;
        t1 = t0 + delta;

        F = zeros(size(t));
        mask0 = (t < t0);
        mask1 = (t >= t0) & (t < t1);
        mask2 = (t >= t1);

        % Smoothstep polynomial (C^2): 6s^5 - 15s^4 + 10s^3
        s = (t(mask1) - t0) / (t1 - t0);
        H = 6*s.^5 - 15*s.^4 + 10*s.^3;

        F(mask1) = Amp * sin(W*t(mask1)) .* H;
        F(mask2) = Amp * sin(W*t(mask2));
        F(mask0) = 0;

    case 2
        % Sigmoid (step-like) ramp starting at t=A
        A = 20;
        k = 0.1;

        F = zeros(size(t));
        mask = (t >= A);
        F(mask) = Amp ./ (1 + exp(-k * (t(mask) - A)));

    case 3
        % C^2 smooth step from 0 to Amp over [A,C]
        A = 20;
        C = 25;

        F = zeros(size(t));
        mask0 = (t < A);
        mask1 = (t >= A) & (t <= C);
        mask2 = (t > C);

        s = (t(mask1) - A) / (C - A);
        H = 6*s.^5 - 15*s.^4 + 10*s.^3;

        F(mask0) = 0;
        F(mask1) = Amp * H;
        F(mask2) = Amp;

    otherwise
        error('MyForce:InvalidFlag','Unknown DATA.Flag value: %g', DATA.Flag);
end

end
