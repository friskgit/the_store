% Length of frame
  s_file = 'bark10_m.wav';
  L = 1024;             
  samples = [10241, 11264];
  [len, fs] = audioread(s_file);
  [x, fs] = audioread(s_file, samples);
  				% sound file length
  sf_length = size(len);
  bin_size = fs/L;
  sf_frames = fix(sf_length / L);
            			% Frequencies in steps of R / fft size
  f = bin_size*(0:(L/2));

  T = 1/fs;             % Sampling period       
  t = (0:L-1)*T;        % Time vector

  				% Window function
  NFFT = 2^nextpow2(L); 
  Y = fft(hanning(L).*x)/L;
          			% Amplitude spectrum in P1
  P2 = abs(Y/L);
  P1 = P2(1:L/2+1);
  P1(2:end-1) = 2*P1(2:end-1);
  				% Calculate filter
  f_bark = arrayfun(@(x) 26.81*x/(1960+x)-0.53, f);
  f_val = 10.^((15.81 + 7.5*(f_bark+0.474) - 17.5*(1 + (f_bark+0.474).^2).^(1/2))/80);
  s_con = f_val' .* P1;

  disp(f_bark(159));
  % Adjust amplitudees
  s_con_sq = sqrt(s_con);
  [mx, imx] = max(s_con_sq);
  				% Plot it
  subplot (2, 1, 1)
  plot(f,P1,"LineWidth",0.5) 
  title("Single-Sided Amplitude Spectrum of X(t)")
  xlabel("f (Hz)")
  ylabel("|amplitude(f)|")

  subplot (2, 1, 2)
  plot(f_bark, s_con_sq,"LineWidth",0.5)
  title("Bark Scale Spectrum of X(t)")
  xlabel("Bark")
  ylabel("|amplitude(bark)|")
%  xlim([1 24])
