s_file = 'bark10_m.wav';
L = 1024;
info = audioinfo(s_file);
[len, fs] = audioread(s_file);
sf_length = info.TotalSamples;
bin_size = fs/L;
sf_frames = fix(sf_length / L);

result = [;];

t = 1:sf_length-1;

f = bin_size*(0:(L/2));

T = 1/fs;             % Sampling period       
t = (0:L-1)*T;        % Time vector
                      % Window function
NFFT = 2^nextpow2(L); 

			      % Loop over the frames of the sound file
for i=1:sf_frames-1
  bin = [];
  start = i*1024;
  stop = start+1023;
  samples = [start, stop];
  [x, fs] = audioread(s_file, samples);

  Y = fft(hanning(L).*x)/L;
                        	% Amplitude spectrum in P1
  P2 = abs(Y/L);
  P1 = P2(1:L/2+1);
  P1(2:end-1) = 2*P1(2:end-1);
                     		% Calculate filter
  f_bark = arrayfun(@(x) 26.81*x/(1960+x)-0.53, f);
  f_val = 10.^((15.81 + 7.5*(f_bark+0.474) - 17.5*(1 + (f_bark+0.474).^2).^(1/2))/80);
          			%    s_con = sqrt(f_val' .* P1);
  s_con = f_val' .* P1;
  s_conq = sqrt(s_con);
  [mx, imx] = max(s_conq);

				% extract bark bands from the fft
  bark_i=0;
  for j=1:512
			 % is the current bark band larger than bark_i
    if f_bark(j) > bark_i;
 				%        disp(f_bark(j));
      bin = [bin, s_conq(j)];
      bark_i = bark_i+1;
    endif
  end
  result = [result;bin];
end
          			% Normalization factor
norm = 1.0/max(max(result));
          			% Normalize
result_norm = norm*result;

test = result_norm(10, :);

plot(test, "LineWidth",0.5)
title("Bark Scale Spectrum of X(t) (normalized)")
xlabel("Bark")
ylabel("|amplitude(bark)|")
xlim([1 24])

file_id = fopen('violin.wav.txt', 'w');
for j=1:40
  fprintf(file_id, '%f ', result_norm(j,:));
fdisp(file_id, ' ');
end
