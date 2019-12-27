function [ mu, v ] = ksmoother( x, alpha, sigma, fc, fs )
%KSMOOTHER Kalman smoother for bandpassed analytical signals
%  x     - analytical signal ([1 T] vector of complex numbers)
%  alpha - variance of measurement noise
%  sigma - variance of signal
%  fc    - center frequency
%  fs    - sampling rate

  w0 = 2 * pi * fc / fs;
  T  = length(x);
  
  % Forward pass
  mu = zeros(1,T);
  p  = zeros(1,T);
  
  mu(1) = x(1);
  p(1)  = alpha;
  r  = exp(1i * w0);
  
  for n = 2:T
      lambda = (p(n-1) + sigma) / (p(n-1) + sigma + alpha);
      mu(n)  = lambda * x(n) + (1 - lambda) * r * mu(n-1);
      p(n)   = alpha * lambda;
  end
  
  % Backward pass
  v = zeros(1,T);
  
  v(T) = p(T);
  ri = exp(-1i * w0);
  
  for n = (T-1):-1:1
      lambda = p(n) / (p(n) + sigma);
      mu(n)  = lambda * ri * mu(n+1) + (1 - lambda) * mu(n);
      v(n)   = p(n) + lambda^2 * (v(n+1) - p(n+1) - sigma);
  end
end

