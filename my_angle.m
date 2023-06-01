function an = my_angle(z)
    if angle(z) < 0
        phi = 2*pi + angle(z);
    else
        phi = angle(z);
    end
    an = phi;
end