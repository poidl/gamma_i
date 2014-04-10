function yout = dj_smuth(yin)

%%%    Usage:    yout = dj_smuth(yin)


ntimes = 20; n = length(yin);

yout = yin;

for itimes = 1:ntimes

  ytmp = 0.25*yout(1:n-2)+0.5*yout(2:n-1)+0.25*yout(3:n);

  yout = [yout(1); ytmp(:); yout(n)];

end

return

