function [out_bits] = demodulator(in_syms, modulation_order)
    if ~(modulation_order==1 || modulation_order==2)
        error('Wrong modulation order!')
    end
    if modulation_order == 1
        out_bits = double(in_syms >= 0);
    else
        temp1 = ([real(in_syms), imag(in_syms)] > 0)*1;
        out_bits = reshape(temp1', [], 1);
    end

end