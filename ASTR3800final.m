clear
%ASTR3800 final project
fifty=importdata("DM_1M.dat");
%linking radius
max(fifty(:,1))
r_l=0.2*(141.3)/(1000000^(1/3));
fifty=fifty.*141.3/100; %put points randomly into the cube

%only need to consider wraparounds for galaxies less than or exactly r_l away from border
fifty_dup=fifty
fifty_dup(fifty_dup>=(141.3-r_l))=fifty_dup(fifty_dup>=(141.3-r_l))-141.3;

%all duplicates
fifty_dup=fifty_dup(fifty_dup(:,1)<0|fifty_dup(:,2)<0|fifty_dup(:,3)<0,:)

%counter for number of halos
pp=0;
ppp=0;
cell_halo=cell(1,500000);

while isempty(fifty)~=1
    
    ppp=0; 
    %reset starting halo
    halo=zeros(length(fifty),3);
    pp=pp+1;
    ppp=ppp+1;
    ppp_original=ppp;

    %first potential start of halo
    halo(1,:)=fifty(1,:);
    fifty(1,:)=[];
    
    %array of distances
    dist_real=(((fifty(:,1)-halo(1,1)).^2)+((fifty(:,2)-halo(1,2)).^2)+((fifty(:,3)-halo(1,3)).^2)).^.5;
  
    dist_duplicate=(((fifty_dup(:,1)-halo(1,1)).^2)+((fifty_dup(:,2)-halo(1,2)).^2)+((fifty_dup(:,3)-halo(1,3)).^2)).^.5;
    
    %indices of distances, add equivalent indices to array
    hal_added=find(dist_real<r_l);
    hal_added_dup=find(dist_duplicate<r_l);
    hal_added_delete=hal_added;
    hal_added_dup_delete=hal_added_dup;
    row_delete=zeros(1,length(hal_added_dup));
    row_count=0;
    %add new galaxies to halos (need to then check all new halos)
    k=[];
    while isempty(hal_added)~=1
        halo(ppp,:)=fifty(hal_added(1),:);
        if fifty(hal_added(1),:)>(141.3-r_l)
            temp_var=hal_added(1);
            fifty(fifty((temp_var),:)>141.3-r_l)=fifty(fifty((temp_var),:)>141.3-r_l)-141.3;
            k=find(fifty(temp_var,1)==fifty_dup(:,1)&fifty(temp_var,2)==fifty_dup(:,2)&fifty(temp_var,3)==fifty_dup(:,3));
        end
        ppp=ppp+1;
        halo(ppp,:)=fifty(hal_added(1),:);
        hal_added(1)=[];
        ppp=ppp+1;
        
    end
    
    while isempty(hal_added_dup)~=1
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
        fifty_dup(fifty_dup(hal_added_dup(1),:)<0)= fifty_dup(fifty_dup(hal_added_dup(1),:)<0)+141.3;
        ppp=ppp+1;
        z=find((fifty(:,1)==fifty_dup((hal_added_dup(1)),1)&fifty(:,2)==fifty_dup((hal_added_dup(1)),2)&fifty(:,3)==fifty_dup((hal_added_dup(1)),3)));
        
        if length(z)>=1
            for ii=1:length(z)
            row_count=row_count+1;
            row_delete(row_count)=z(ii)
            end
            
        else
            row_delete=[];
        end
        row_count=row_count+1;
        halo(ppp,:)=fifty_dup(hal_added_dup(1),:);
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
 
 if isempty(row_delete)==0
     fifty(row_delete(:),:)=[];
 end
    %halos starting from 1,1
    
    numb_added=ppp-ppp_original;
 
    %need to eliminate duplicates here in code maybe
 while numb_added~=0

    
  
   
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
    %add new galaxies to halos (need to then check all new halos)
    
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
 
 ind=find(sum(halo,2)==0);
 halo(ind,:)=[];
halo(halo<0)=halo(halo<0)+141.3;
halo=unique(halo, 'rows');
    %at this point, all galaxies both real and duplicate a distance of r_l
    %away from og galaxy recorded, and all galaxies that distance from

    %first set of galaxies recorded
cell_halo{pp}=halo;
end
cell_halo=cell_halo(~cellfun('isempty',cell_halo));





isempty(fifty)
%could use case, where each cell array is a halo, and each entry in cell
%array is xyz coordinate, then for each halo, see from
cell_halo