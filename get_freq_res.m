function [w, freq_res_og, freq_res_test] = get_freq_res(og_zeros, test_zeros)

w = linspace(-pi,pi,1024);
z = exp(1i*w);

freq_res_og = ones(1,length(w));
freq_res_test = ones(1,length(w));

for index = 1:length(og_zeros)
    freq_res_og = freq_res_og.*(z - og_zeros(index))./z;
    freq_res_test = freq_res_test.*(z - test_zeros(index))./z;
end
end