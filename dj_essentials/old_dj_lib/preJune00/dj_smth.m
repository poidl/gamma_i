function yout = dj_smth(yin)

%%%    Usage:    yout = dj_smth(yin)
%%%

ntimes = 5; n = length(yin);

yout = yin;

for itimes = 1:ntimes

  ytmp = 0.1*yout(1:n-4)+0.2*yout(2:n-3)+0.4*yout(3:n-2)+ ...
            0.2*yout(4:n-1)+0.1*yout(5:n);

  yout = [yout(1); yout(2); ytmp; yout(n-1); yout(n)];

end

return

