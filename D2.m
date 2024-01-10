function df = D2(r,rh,s)
df = 1./sqrt(r.^2+s.^2)-(rh.^2)./(r.^2+s.^2).^(1.5);
% df = (1+2.*sqrt(r.^2))+2.*sqrt(r.^2);
% df = (-2*(s^2)+4*(s^4).*r.^2).*exp(-(s.*r).^2);
end