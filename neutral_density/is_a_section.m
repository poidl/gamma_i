function answer = is_a_section(s,t)

dims = size(s); answer = 0;

if dims(3)==3
    ss = permute(s,[3,1,2]); tt = permute(t,[3,1,2]);
    zs = diff(ss); zt = diff(tt); 
    answer = 1-(nanmax(abs(zs(:)))+nanmax(abs(zt(:))));
end

if answer~=1, answer = 0; end

return