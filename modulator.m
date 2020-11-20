function out_syms = modulator(in_bits, modulation_order)

if (modulation_order~=1 && modulation_order~=2)
    error('Modulation_order wrong!')
end

if modulation_order == 1 % BPSK
    out_syms = 1-2*in_bits; % 0 -> 1; 1 -> -1 
else
   temp2 = 2*reshape(in_bits, 2, []).'-1;
   out_syms = (1/sqrt(2))*(temp2(:, 1)+1i*temp2(:, 2));
end

end