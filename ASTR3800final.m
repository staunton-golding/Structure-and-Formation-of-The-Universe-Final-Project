clear
%ASTR3800 final project
%fifty and all variants are the variable names regardless of size of
%dataset as the first dataset considerd was the 50k one
fifty=importdata("DM_200k.dat");
%data is normalized to max values = 100, need to renormalize to max =
%141.3
fifty=fifty.*141.3/100
%find linking radius
r_l=0.2*(141.3)/(200000^(1/3));

%duplicate array of all DM
fifty_dup=fifty;
%all DM less than linking radius away from border needs to be duplicated to
%account for wraparounds
fifty_dup(fifty_dup>=(141.3-r_l))=fifty_dup(fifty_dup>=(141.3-r_l))-141.3;

%only needed duplicates are kept
fifty_dup=fifty_dup(fifty_dup(:,1)<0|fifty_dup(:,2)<0|fifty_dup(:,3)<0,:)

%counter for number of halos
pp=0;

%cell to hold array of all halos (for 1M, less than 500,000 halos,
%arbitrary value based on literature)
cell_halo=cell(1,500000);

%while DM matter still unassigned to halos, halos still to be formed
while isempty(fifty)~=1
    %counter
    ppp=0; 
    %to preallocate for speed, no halo is going to be bigger than all
    %remaining DMs
    halo=zeros(length(fifty),3);
    %counter increase (for cell_halo indexing)
    pp=pp+1;
    %to start halo entries at 1
    ppp=ppp+1;
    %original length of halo (later used to track increase in halo size)
    ppp_original=ppp;
    %first potential start of halo
    halo(1,:)=fifty(1,:);
    %eliminate DM added to halo
    fifty(1,:)=[];
    
    %array of distances from DM that starts the perspective halo
    dist_real=(((fifty(:,1)-halo(1,1)).^2)+((fifty(:,2)-halo(1,2)).^2)+((fifty(:,3)-halo(1,3)).^2)).^.5;
  
    %array of distances from DM when accounting for potential wraparounds
    dist_duplicate=(((fifty_dup(:,1)-halo(1,1)).^2)+((fifty_dup(:,2)-halo(1,2)).^2)+((fifty_dup(:,3)-halo(1,3)).^2)).^.5;
    
    %indices DM's within linking radius from original DM
    hal_added=find(dist_real<r_l);
    %indices of duplicated DM's within linking radius from original DM
    hal_added_dup=find(dist_duplicate<r_l);
    %indices of DM's to be deleted after they are added to the halo
    hal_added_delete=hal_added;
    %indices of duplicates DM's to be deleted after they are added to the
    %halo
    hal_added_dup_delete=hal_added_dup;
    
    %DM's that potentially translate to (+) counterparts from duplicate arrays need to be kept
    %track of for deletion within original array
    row_delete=zeros(1,length(hal_added_dup));
    %number of row's from duplicate to be deleted (starts at zero)
    row_count=0;
    %indices for row_count, starts as empty until a DM added from original
    %array has a duplicated counterpart
    k=[];

    %adds DM to original halo
    while isempty(hal_added)~=1
        %added DM corresponds to first entry of hal_added- an array holding
        %indices of DM to be added to halo
        halo(ppp,:)=fifty(hal_added(1),:);
        %check to see if DM added has duplicated counterpart
        if fifty(hal_added(1),:)>(141.3-r_l)
            temp_var=hal_added(1);
            %set the row in original array to = duplicated counterpart
            fifty(fifty((temp_var),:)>141.3-r_l)=fifty(fifty((temp_var),:)>141.3-r_l)-141.3;
            %find index of duplicated counterpart in duplicated array
            k=find(fifty(temp_var,1)==fifty_dup(:,1)&fifty(temp_var,2)==fifty_dup(:,2)&fifty(temp_var,3)==fifty_dup(:,3));
        end
        %counter
        ppp=ppp+1;
        %added the negative transfrom duplciated counterpart to halo (will
        %later be removed, accounts for wraparounds)
        halo(ppp,:)=fifty(hal_added(1),:);
        %delete first entry in array holding indices of added DM
        hal_added(1)=[];
        %counter
        ppp=ppp+1;  
    end
 
    
    %add duplicate halos
    while isempty(hal_added_dup)~=1
        %add duplicate halo w/ index from hal_added_dup(licate)
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
        %find the (+) counterpart to duplicated DM
        fifty_dup(fifty_dup(hal_added_dup(1),:)<0)= fifty_dup(fifty_dup(hal_added_dup(1),:)<0)+141.3;
        %counter
        ppp=ppp+1;
        %find index of (+) original equal to duplicate
        z=find((fifty(:,1)==fifty_dup((hal_added_dup(1)),1)&fifty(:,2)==fifty_dup((hal_added_dup(1)),2)&fifty(:,3)==fifty_dup((hal_added_dup(1)),3)));
        %if the original counterpart exists, and if its >1 (unlikely but for robustness)
        if length(z)>=1
            %for all z's
            for ii=1:length(z)
            %counter
            row_count=row_count+1;
            %add index of row to be deleted from original array
            row_delete(row_count)=z(ii)
            end
        else
            %if z=0 or [], there will be no row deletions from original
            %array
            row_delete=[];
        end
        %counter (might be duplicated rows later)
        row_count=row_count+1;
        %add original counterpart to duplicated DM
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
        %delete index of first DM added from duplicate array
        hal_added_dup(1)=[];
        %counter
        ppp=ppp+1;
    end
 %delete all added DMS from original array
 if isempty(hal_added_delete)==0    
    fifty(hal_added_delete(:),:)=[];
 end

 %delete all DM added to halo from duplicated array
 if isempty(hal_added_dup_delete)==0 && isempty(k)==0
     %both the duplicated and duplicate counterparts added
     hal_added_dup_delete=[hal_added_dup_delete, k];
     %all nonunique entries deleted (unlikely)
     hal_added_dup_delete=unique(hal_added_dup_delete);
     %all added DM's from duplicate or duplicate counterpart deleted
     fifty_dup(hal_added_dup_delete(:),:)=[];
 elseif isempty(hal_added_dup_delete)==0 && isempty(k)==1
     fifty_dup(hal_added_dup_delete(:),:)=[];
 else 
     fifty_dup(k(:),:)=[];
 end

 %all original counterparts to duplicated DM additions deleted
 if isempty(row_delete)==0
     fifty(row_delete(:),:)=[];
 end
   
 %number of added DMs to halo
 numb_added=ppp-ppp_original;
 
 %need to check each new addition for DM less than linking radius away
 while numb_added~=0
    %if only one DM is added, repeat lines 35 through 130
  if numb_added==1
    ppp_or=ppp;
        
    dist_real=(((fifty(:,1)-halo(ppp,1)).^2)+((fifty(:,2)-halo(ppp,2)).^2)+((fifty(:,3)-halo(ppp,3)).^2)).^.5;
    dist_duplicate=(((fifty_dup(:,1)-halo(ppp,1)).^2)+((fifty_dup(:,2)-halo(ppp,2)).^2)+((fifty_dup(:,3)-halo(ppp,3)).^2)).^.5;
    %then need indices of distances, add equivalent indices to array
    hal_added=find(dist_real<r_l);
    hal_added_dup=find(dist_duplicate<r_l);
    hal_added_delete=hal_added;
    hal_added_dup_delete=hal_added_dup;
    row_delete=zeros(1,length(hal_added_dup));
    row_count=0;
    k=[];
    
    while isempty(hal_added)~=1
        halo(ppp,:)=fifty(hal_added(1),:);
        if fifty(hal_added(1),:)>(141.3-r_l)
            temp_var=hal_added(1);
            fifty(fifty((temp_var),:)>141.3-r_l)=fifty(fifty((temp_var),:)>141.3-r_l)-141.3;
            k=find(fifty(temp_var,1)==fifty_dup(:,1)&fifty(temp_var,2)==fifty_ddup(:,2)&fifty(temp_var,3)==fifty_dup(:,3));
        end
        ppp=ppp+1;
        halo(ppp,:)=fifty(hal_added(1),:);
        ppp=ppp+1;
        hal_added(1)=[];
    end

    while isempty(hal_added_dup)~=1
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
        ppp=ppp+1;
        fifty_dup(fifty_dup(hal_added_dup(1),:)<0)= fifty_dup(fifty_dup(hal_added_dup(1),:)<0)+141.3;
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
        z=find((fifty(:,1)==fifty_dup((hal_added_dup(1)),1)&fifty(:,2)==fifty_dup((hal_added_dup(1)),2)&fifty(:,3)==fifty_dup((hal_added_dup(1)),3)));
        if length(z)>=1
            for ii=1:length(z)
            row_count=row_count+1;
            row_delete(row_count)=z(ii)
            end
        else
            row_delete=[];
        end
        hal_added_dup(1)=[];
        ppp=ppp+1;
    end

    if isempty(hal_added_delete)==0    
       fifty(hal_added_delete(:),:)=[];
    end
    if isempty(hal_added_dup_delete)==0 && isempty(k)==0
       hal_added_dup_delete=[hal_added_dup_delete, k];
       hal_added_dup_delete=unique(hal_added_dup_delete);
       fifty_dup(hal_added_dup_delete(:),:)=[];
    elseif isempty(hal_added_dup_delete)==0 && isempty(k)==1
        fifty_dup(hal_added_dup_delete(:),:)=[];
    else 
       fifty_dup(k(:),:)=[];
    end 
    if isempty(row_delete)~=1
       fifty(row_delete(:),:)=[];  
    end
    numb_added=ppp-ppp_or;
 else
    ppp_or=ppp;
    jay=ppp;
     
%if more than one added, sets each added DM as the base DM to check distances for new halo additions.
%This repeats for every new addition added from above. Once each of the new entries are accounted for, if even more entries are added, this proccess is repeated until all DM's added to halo are checked
%for DM's less than a linking radius away. Process terminates when all new
%entries are checked, and no new DM's are added to the halo.

    %same steps as lines 35-130 and lines 173 through 217
    for kk=(jay-numb_added):jay         
        dist_real=(((fifty(:,1)-halo(kk,1)).^2)+((fifty(:,2)-halo(kk,2)).^2)+((fifty(:,3)-halo(kk,3)).^2)).^.5;
        dist_duplicate=(((fifty_dup(:,1)-halo(kk,1)).^2)+((fifty_dup(:,2)-halo(kk,2)).^2)+((fifty_dup(:,3)-halo(kk,3)).^2)).^.5;
        hal_added_new=find(dist_real<r_l);
        hal_added_dup_new=find(dist_duplicate<r_l);
        hal_added_delete_new=hal_added_new;
        hal_added_dup_delete=hal_added_dup_new;
        row_delete=zeros(1,length(hal_added_dup));
        row_count=0;
        k=[];
        while isempty(hal_added_new)~=1
           halo(ppp,:)=fifty(hal_added_new(1),:);
           if fifty(hal_added_new(1),:)>(141.3-r_l)
               temp_var=hal_added_new(1);
               fifty(fifty((temp_var),:)>141.3-r_l)=fifty(fifty((temp_var),:)>141.3-r_l)-141.3;
               k=find(fifty(temp_var,1)==fifty_dup(:,1)&fifty(temp_var,2)==fifty_dup(:,2)&fifty(temp_var,3)==fifty_dup(:,3));
           end
           ppp=ppp+1;
           halo(ppp,:)=fifty(hal_added_new(1),:);
           ppp=ppp+1;
           hal_added_new (1)=[];
        end

        while isempty(hal_added_dup_new)~=1    
           halo(ppp,:)=fifty_dup(hal_added_dup_new(1),:);
           ppp=ppp+1;
           fifty_dup(fifty_dup(hal_added_dup_new(1),:)<0)= fifty_dup(fifty_dup(hal_added_dup_new(1),:)<0)+141.3;       
           z=find((fifty(:,1)==fifty_dup((hal_added_dup_new(1)),1)&fifty(:,2)==fifty_dup((hal_added_dup_new(1)),2)&fifty(:,3)==fifty_dup((hal_added_dup_new(1)),3)));
            
           if length(z)>=1
               for ii=1:length(z)
               row_count=row_count+1;
               row_delete(row_count)=z(ii);
               end
           else
               row_delete=[];
           end
           
           halo(ppp,:)=fifty_dup(hal_added_dup_new(1),:);
           ppp=ppp+1;
           hal_added_dup_new(1)=[];
        end
         
        if isempty(hal_added_delete_new)==0    
            fifty(hal_added_delete_new(:),:)=[];
        end

        if isempty(hal_added_dup_delete)==0 && isempty(k)==0
            hal_added_dup_delete=[hal_added_dup_delete, k];
            hal_added_dup_delete=unique(hal_added_dup_delete);
            fifty_dup(hal_added_dup_delete(:),:)=[];
        elseif isempty(hal_added_dup_delete)==0 && isempty(k)==1
            fifty_dup(hal_added_dup_delete(:),:)=[];
        else 
            fifty_dup(k(:),:)=[];
        end

        if isempty(row_delete)==0
           fifty(row_delete(:),:)=[];
        end
        numb_added=ppp-ppp_or;
    end 
  end
end

%finds all zero entries into the halo, and eliminate them
ind=find(sum(halo,2)==0);
halo(ind,:)=[];
%all duplicates are made original again (+)
halo(halo<0)=halo(halo<0)+141.3;
%delete all repeated rows
halo=unique(halo, 'rows');
    
%enters halo array into cell that holds each individual halo
cell_halo{pp}=halo;
end

%eliminates all 'empty' entries into the preallocated cell_halo cell array
cell_halo=cell_halo(~cellfun('isempty',cell_halo));
%counter
pp=0;
%checks for all single DM halos from cell_halo

for jj=1:length(cell_halo)
    [row col]=size(cell_halo{jj});
    if row==1
        pp=pp+1;
        cell_2_delete(pp)=jj;
    end
end

%makes all single DM halos empty
for ii=1:length(cell_2_delete)
    cell_halo{cell_2_delete(ii)}=[];
end

%deletes all empty (once single DM) halos
cell_halo=cell_halo(~cellfun('isempty',cell_halo));

%to hold the rms Radius for all halos
rad=zeros(1,length(cell_halo));
%to hold total masses of all halos
total_mass=zeros(1,length(cell_halo));
%to hold center of masses for all halos
center_om=zeros(3, length(cell_halo));

%to find rad, total_mass and center_om, plus account for wrapping 
for jj=1:length(cell_halo)
    %halo considered in each cycle of loop
    temp_hal=cell_halo{jj};
    %size of each halo
    [row col]=size(temp_hal);
    %total mass of each halo
    total_mass(jj)=row*1.4*10^10;
    %x correction for center of mass and radius
    %distance min x coordiante is away from 0
    x_min=min(temp_hal(:,1));
    x_max=max(temp_hal(:,1));
    %distance max x coordinate is away from 141.3
    dist_bound=141.3-x_max;
    %range of x coordinates
    x_range=x_max-x_min;
    
    %if the range of x-coordinates is greater than 141.3/2, then there will
    %be wrapping
    if x_range >141.3/2
        %if data is closer to 0 than to 141.3, data farther than 141.3/2 from x_min needs to be wrapped
        %around 141.3. Opposite if data is closer to 141.3 than 0
        if dist_bound > x_min
            temp_hal(temp_hal(:,1)<(dist_bound-141.3/2))=temp_hal(temp_hal(:,1)<(dist_bound-141.3/2))+141.3;
        else
            temp_hal(temp_hal(:,1)>(x_min+141.3/2))=temp_hal(temp_hal(:,1)>(x_min+141.3/2))-141.3;
        end
    end
    
    
    %y correction (repeat lines 346 to 363 except for y)
    y_min=min(temp_hal(:,2));
    y_max=max(temp_hal(:,2));
    dist_bound=141.3-y_max;
    y_range=y_max-y_min;
    
    if y_range >141.3/2
        if dist_bound > y_min
            temp_hal(temp_hal(:,2)<(dist_bound-141.3/2))=temp_hal(temp_hal(:,2)<(dist_bound-141.3/2))+141.3;
        else
            temp_hal(temp_hal(:,2)>(y_min+141.3/2))=temp_hal(temp_hal(:,2)>(y_min+141.3/2))-141.3;
        end
    end
    %z correction (repeat liens 346 to 363 except for z)
    z_min=min(temp_hal(:,1));
    z_max=max(temp_hal(:,1));
    dist_bound=141.3-z_max;
    z_range=z_max-z_min;
    
    if z_range >141.3/2
        if dist_bound > z_min
            temp_hal(temp_hal(:,3)<(dist_bound-141.3/2))=temp_hal(temp_hal(:,3)<(dist_bound-141.3/2))+141.3;
        else
            temp_hal(temp_hal(:,3)>(z_min+141.3/2))=temp_hal(temp_hal(:,3)>(z_min+141.3/2))-141.3;
        end
    end
    %as all masses same, average coordinates = center of mass coordinates
    x_center=sum(temp_hal(:,1))/row;
    y_center=sum(temp_hal(:,2))/row;
    z_center=sum(temp_hal(:,3))/row;
    %array to hold radii of all DM from center of mass fro each halo
    radii=zeros(1,row);
    for nn=1:row
        %find radii
        radii(nn)=((temp_hal(nn,1)-x_center)^2+(temp_hal(nn,2)-y_center)^2+(temp_hal(nn,3)-z_center)^2)^.5;
    end
    %square all radii radius
    radii=radii.^2;
    %sum radii
    rad(jj)=sum(radii);
    %take squareroot of sum to find RMS value
    rad(jj)=rad(jj)^.5;
end

%to find number density of halos greater than a given mass
n_M=zeros(1,1000);
%mass thresholds considered
mass=(0:max(total_mass)/1000:max(total_mass));
%counter
pp=0;
for jj=0:max(total_mass)/1000:max(total_mass)
    %counter
    pp=pp+1;
    %find all halos greater than a given mass
    x=find(total_mass(:)>jj);
    %find number density of halos greater than a given mass
    n_M(pp)=length(x)/(141.3^3);
end

%scatter plots
scatter(log10(mass),log10(n_M))
title('log10(Mass) vs log10(n(>M))')
xlabel('log10(Mass)')
ylabel('log10(n(>M))')

scatter(log10(total_mass),log10(rad))
title('log10(Mass) vs log10(RMS radius)')
xlabel('log10(Mass)')
ylabel('log10(RMS radius)')
