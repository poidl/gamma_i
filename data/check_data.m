
%%      deletes casts with internal missing data and less than three levels
                           

inds = find(isfinite(s(1,:))); [nz,nxy] = size(s);
ss = s(:,inds); z = flipud(ss); [nzz,nxx] = size(z);
zz = nan*ones(1,nxx);
for i = 1:nxx
    indss = find(isfinite(z(:,i))); zz(i) = nz-indss(1);
end


%   (i)

inds_delete = find(n(inds)<zz);

deleting_casts = inds_delete, salt = s(:,inds(inds_delete));

s(:,inds(inds_delete)) = nan; t(:,inds(inds_delete)) = nan; g(:,inds(inds_delete)) = nan;

ocean(inds(inds_delete)) = nan; n(inds(inds_delete)) = nan;


%   (ii)

inds_delete = find(n<3); 

number_of_shallow_casts = length(inds_delete), %salt = s(:,inds_delete)

s(:,inds_delete) = nan; t(:,inds_delete) = nan; g(:,inds_delete) = nan;

ocean(inds_delete) = nan; n(inds_delete) = nan;



