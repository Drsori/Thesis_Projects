function df = D1(r,rh,s)
df = rh./sqrt(r.^2+s.^2);
% df = rh.*(1+2.*sqrt(r.^2));
% df = (-2*(s^2).*rh).*exp(-(s.*r).^2);
end