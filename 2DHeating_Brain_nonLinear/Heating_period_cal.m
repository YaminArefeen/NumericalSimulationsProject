function [ Period_acc ] = Heating_period_cal( Temperature,t,temp_thresh,mask,Time_limit_tumor, Time_limit_health)

Period_acc=zeros(size(Temperature,1),1);
for i=2:size(t,1)
    tmp=(Temperature(:,i)+Temperature(:,i-1))/2;
    Period_acc=Period_acc+(tmp>temp_thresh).*(t(i)-t(i-1));
end


whole_periodmap=zeros(size(mask));

whole_periodmap(mask(:)>0)=Period_acc;
killed_tumor=(whole_periodmap>=Time_limit_tumor).*(mask==1);
nonkilled_tumor=(whole_periodmap<Time_limit_tumor).*(mask==1);

hurt_health=(whole_periodmap>=Time_limit_health).*(mask>1);
nonhurt_health=(whole_periodmap<Time_limit_health).*(mask>1);

figure;
im_color_show=zeros([size(mask),3]);
im_color_show(:,:,3)=im_color_show(:,:,3)+logical(nonhurt_health);  % non-hurt health blue
im_color_show(:,:,1)=im_color_show(:,:,1)+logical(hurt_health);
im_color_show(:,:,3)=im_color_show(:,:,3)+logical(hurt_health);      % hurt_health purple
im_color_show(:,:,2)=im_color_show(:,:,2)+logical(killed_tumor);     % killed_tumor green
im_color_show(:,:,1)=im_color_show(:,:,1)+logical(nonkilled_tumor);  % non-killed tumor red
imshow(im_color_show);

end

