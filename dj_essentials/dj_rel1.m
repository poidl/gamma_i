function v = dj_rel1(tls1,wow1)

v = 260*(wow1-11.98) - 600*(tls1-5.05);

v = round(100*v)/100;

return
