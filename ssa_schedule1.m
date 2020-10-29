function [foodfitness,foodposition,convergence_curve]=ssa_schedule1(N,Max_iter,lb,ub,t)
% N : jumlah salp
% Max_iter :Iterasi maksimum
% lb : Lower Bound for int. of the swarm 
% ub : Upper Bound for init. of the swarm
% t = matriks t (tetap diketik t. sebab matriks telah diinput diawal)

ct=cputime;
[job , mesin]= size (t); 
% size(t) : ukuran matrik t
dim=job;% dimensi masalah = jumlah job
if size(ub,1)=1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
    % ones(dim,1)*ub : matrik dengan anggota ub berukuran (dim,1)
end

convergence_curve = zeros(1,Max_iter);
% elemen matrik bernilai 0 dengan ukuran (1,iterasi maksimum)

%Initialize the positions of salps
salppositions=initialization(N,dim,ub,lb);


foodposition=zeros(1,dim);
%posisi makanan : matrik bernilai 0 dengan ukuran(1,job/dim)
foodfitness=inf; % inf : tidak terbatas


%calculate the fitness of initial salps

for i=1:size(salppositions,1)
    salpfitness(1,i)=flowshop(salppositions(i,:),t);
end

[sorted_salps_fitness,sorted_indexes]=sort(salpfitness);

for newindex=1:N
    sorted_salps(newindex,:)=salppositions(sorted_indexes(newindex),:);
end

foodposition=sorted_salps(1,:);
foodfitness=sorted_salps_fitness(1);

%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while l<Max_iter+1
    
    c1 = 2*exp(-(4*l/Max_iter)^2); % Eq. (3.2) in the paper
    
    for i=1:size(salppositions,1)
        
        salppositions= salppositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    salppositions(j,i)=foodposition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    salppositions(j,i)=foodposition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=salppositions(:,i-1);
            point2=salppositions(:,i);
            
            salppositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        salppositions= salppositions';
    end
    
    for i=1:size(salppositions,1)
        
        Tp=salppositions(i,:)>ub';Tm=salppositions(i,:)<lb';salppositions(i,:)=(salppositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        
        salpfitness(1,i)=flowshop(salppositions(i,:),t);
        
        if salpfitness(1,i)<foodfitness
            foodposition=salppositions(i,:);
            foodfitness=salpfitness(1,i);
            
        end
    end
    
    convergence_curve(l)=foodfitness;
    l = l + 1;
end
waktukomputasi=cputime-ct;



