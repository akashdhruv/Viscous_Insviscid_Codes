function [Hstari] = calc_hstar(Hi,Na)

in = length(Hi);
Hstari = zeros(in,1);

for i = 1:in
    if i >Na
        if lt(Hi(i),3.5)
            Hstari(i) = 1.50+0.025*(3.5*Hi(i))^3+0.001*(3.5-Hi(i))^5;
        else
            Hstari(i) = 1.50+0.015*((Hi(i)-3.5)^2)/Hi(i);
        end
    else
        if lt(Hi(i),4)
            Hstari(i) = 1.515+0.076*((Hi(i)-4)^2)/Hi(i);
        else
            Hstari(i) = 1.515+0.040*((Hi(i)-4)^2)/Hi(i);
        end
    end
end

end