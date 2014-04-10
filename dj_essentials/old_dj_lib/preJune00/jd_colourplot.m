% COLOURPLOT  x-y plot of symbols, whose colour is determined by z.
% INPUT:  xx    x coords
%         yy    y coords
%         zz    z coordinate which determines colour
% Optional:       (can use empty [] if require some but not others) 
%         symb   string defining a normal "plot" symbol
%         cpmin  Min value for colour map
%         cprng  Range of colour map (so caxis will be [cpmin cpmin+cprng])
%         msz    Marker size for symbols
%
% Note: Line marker symbols are catered for, but may consume much more time
%       and memory.
%
% Must have previously applied the required colour map.
%
% Jeff Dunn  8/11/96  15/1/98  19/6/98
%
% USAGE: colourplot(xx,yy,zz,symb,cpmin,cprng,msz)

function colourplot(xx,yy,zz,symb,cpmin,cprng,msz)

if nargin==0
  disp('colourplot(xx,yy,zz,{symb,cpmin,cprng,msz})');
  return
end

hold on

rej = find(isnan(zz));
if ~isempty(rej)
  xx(rej) = [];
  yy(rej) = [];
  zz(rej) = [];
end

if min(size(xx)) > 1
  xx = xx(:);
  yy = yy(:);
  zz = zz(:);
end

if nargin < 4 | isempty(symb)
  symb = '+';
end
if nargin < 5 | isempty(cpmin)
  cpmin = min(zz);
end
if nargin < 6 | isempty(cprng)
  cprng = max(zz) - cpmin;
end
if nargin < 7
  msz = [];
end

cm = colormap;
lcm = size(cm,1);



if strcmp(symb(1),'-') | strcmp(symb,':')
  
  % Plot data in original order, esp. for drawing continuous lines

  limit = find(zz<cpmin);
  if ~isempty(limit)
    zz(limit) = cpmin*ones(size(limit));
    disp([num2str(length(limit)) ' data increased to CPMIN'])
  end

  limit = find(zz>(cpmin+cprng));
  if ~isempty(limit)
    zz(limit) = (cpmin+cprng)*ones(size(limit));
    disp([num2str(length(limit)) ' data decreased to CPMIN+CPRNG'])
  end

  cc = cm(1+ floor((lcm-1)*(zz-cpmin)./cprng),:);

  for ii=2:length(xx)
    plot(xx([ii-1 ii]),yy([ii-1 ii]),symb,'Color',cc(ii,:));
  end

else
  
  % Plot in groups of like value, to reduce time and memory consumption.
  
  idx = 1:length(zz);

  if ~isempty(msz)
  
    for ii=1:lcm-1
      kk = find(lcm*(zz(idx)-cpmin)./cprng < ii);
      plot(xx(idx(kk)),yy(idx(kk)),symb,'Color',cm(ii,:),'MarkerSize',msz);
      idx(kk) = [];
    end
    plot(xx(idx),yy(idx),symb,'Color',cm(lcm,:),'MarkerSize',msz);

  else

    for ii=1:lcm-1
      kk = find(lcm*(zz(idx)-cpmin)./cprng < ii);
      plot(xx(idx(kk)),yy(idx(kk)),symb,'Color',cm(ii,:));
      idx(kk) = [];
    end
    plot(xx(idx),yy(idx),symb,'Color',cm(lcm,:));

  end
end

% And now, a terrible slight-of-hand to get colorbar to work properly -
% the last child of the axis must be a patch...
    
caxis([cpmin cpmin+cprng])
aa = axis;
h = patch(aa(1),aa(3),[1 1 1]);
set(h,'EdgeColor','none');

% --------------- End of colourplot.m ----------------------
