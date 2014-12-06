function [new_psi] = apply_BSH_2D_FFT(Vpsi, K2, mu)
Vk = reshape(Vpsi, size(K2));
Vk = fftshift(fft2(Vk));
Vk = ifftshift(Vk./(K2 + mu.^2));
%new_psi = -2.0*real(ifft2(Vk));
new_psi = -2.0*(ifft2(Vk));
new_psi = reshape(new_psi, size(Vpsi));
end

