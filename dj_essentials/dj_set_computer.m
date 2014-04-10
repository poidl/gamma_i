function my_computer = set_computer

[status,cname] = dos('set COMPUTERNAME');

cname = cname(14:23);

if strcmp('C000290-HF',cname)>0
    my_computer = 'desktop';
elseif strcmp('C000068-HF',cname)>0
    my_computer = 'laptop';
else
    my_computer = 'not yet implemented'
end

%the_computer = [cname,'   ', my_computer];

return