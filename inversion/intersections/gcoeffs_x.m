  
gc1a = -(dsige*dsigl*dsigu*(longs(j,i) - longs(j,1 + i))*(pp(-1 + k) - pp(1 + k)));
  
gc1b = dsigl*dsigu*dsigw*(longs(j,-1 + i) - longs(j,i))*(pp(-1 + k) - pp(1 + k));
  
gc1c = -(dsige*dsigw*(longs(j,-1 + i) - longs(j,1 + i))*(dsigu*pp(-1 + k) - ...
                                           dsigl*pp(k) - dsigu*pp(k) + dsigl*pp(1 + k)));
  
gc1(kg) =  gc1a + gc1b + gc1c;
     
gc2(kg) = -(dsige*dsigl*dsigw*(longs(j,-1 + i) - longs(j,1 + i))*(pp(k) - pp(1 + k)));
     
gc3(kg) = dsige*dsigu*dsigw*(longs(j,-1 + i) - longs(j,1 + i))*(pp(-1 + k) - pp(k));
     
gc4(kg) =  -(dsigl*dsigu*dsigw*(longs(j,-1 + i) - longs(j,i))*(pp(-1 + k) - pp(1 + k)));
     
gc5(kg) =  dsige*dsigl*dsigu*(longs(j,i) - longs(j,1 + i))*(pp(-1 + k) - pp(1 + k));


% out = [dsigl,dsigu,dsige,dsigw]
% 
% out = [k,j,i,longs(j,-1 + i),longs(j,1 + i),pp(-1 + k),pp(k)]
% 
% out = [kg,gc1(kg),gc2(kg),gc3(kg),gc4(kg),gc5(kg)], dj_pause(0)