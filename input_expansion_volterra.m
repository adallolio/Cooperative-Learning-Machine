function exp_input = input_expansion_volterra(Mi, exord, F)

if exord >= 1
    h1 = zeros(1,Mi)
    for i = 1:Mi
        h1(i) = F.xBuff(1+i);
    end
end
if exord >= 2
    h2 = zeros(1,Mi^2)
    for i = 1:Mi
        for j = 1:Mi
            h2(i) = F.xBuff(1+i)*F.xBuff(1+j);
        end
    end
end
if exord >= 3
    h3 = zeros(1,Mi^3)
    for i = 1:Mi
        for j = 1:Mi
            for k = 1:Mi
                h3(i) = F.xBuff(1+i)*F.xBuff(1+j)*F.xBuff(1+k);
            end
        end
    end
end

exp_input = horzcat(h1,h2,h3);

end

