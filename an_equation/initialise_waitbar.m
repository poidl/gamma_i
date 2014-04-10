handle_waitbar =  waitbar(0,['generated ...   0%'],'Position',[561 479 270 47]); 

handle_waitbar_c1 = get(handle_waitbar,'Children');

handles_waitbar_c2 = get(handle_waitbar_c1,'Children');

set(handles_waitbar_c2(2),'FaceColor',[0 0 1],'EdgeColor',[0 0 0.5])