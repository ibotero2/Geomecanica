function  [v_Forces, NormalF_p1, TangF_p1  ] = contForces2(cc_part_i, U_part_i, V_part_i, MatPart, UPart, VPart, NormalF, TangF, ks, kn, ro, dt, partIDs)

[mm nn]= size(MatPart);

ll = length(NormalF);
NormalF_p1 = zeros(ll,1);
TangF_p1   = zeros(ll,1);

v_Forces = zeros(3,mm);

if (nnz(UPart)==0 && nnz(U_part_i)==0)
    NormalF_p1 = NormalF;
    TangF_p1  = TangF;
    return;
end


for i=1:mm
   
    cc_part_j = MatPart(i,:);
    U_part_j  = UPart(i,:);
    V_part_j  = VPart(i,:);
    
    
    if ( (cc_part_j(1)==cc_part_i(1)) && (cc_part_j(2)==cc_part_i(2)) )
        v_Forces(:,i) = 0;
    else
       % get n-s local trans matriz
       d2= (cc_part_j-cc_part_i)*(cc_part_j-cc_part_i)';
       
       Fnt = NormalF(partIDs(i));
       Fst = TangF(partIDs(i));
       
       if sqrt(d2)<=2*ro
           e= (cc_part_j-cc_part_i)'/sqrt(d2);
           
           vn = e(1)* (V_part_i(1) - V_part_j(1)) + e(2)* (V_part_i(2) - V_part_j(2));
           vs = e(2)* (V_part_i(1) - V_part_j(1)) - e(1)* (V_part_i(2) - V_part_j(2)) - ro*(V_part_i(3) + V_part_j(3))/1000;
           
           dfn  = kn * vn * dt;
           dfs  = ks * vs * dt;
           
           Fn_ij = Fnt + dfn;
           Fs_ij = Fst + dfs;
           
           
           if (Fn_ij<0)
               Fn_ij=0;
               Fs_ij=0;
           end
           
%            Fs_ijmax = Fn_ij*tan(30*pi/180) + 50000*1e-7*1*ro*(30*pi/180);
           Fs_ijmax = Fn_ij*tan(30*pi/180);
           
           if abs(Fs_ij)>Fs_ijmax
               Fs_ij = abs(Fs_ijmax)*sign(Fs_ij);
           end
           
           
           NormalF_p1(partIDs(i)) = Fn_ij;
           TangF_p1(partIDs(i))  =  Fs_ij;
           
           
                                 
% % %            floc(1) = -abs(floc(1));
% % %            floc(2) = 0;
% % %            
% % %                
% % %                MC = abs(floc(1))*tan(pi/6);
% % %                if abs(floc(2))>MC
% % %                    floc(2) = MC*sign(floc(2));
% % %                end
               
           beta_i = -[e(1) e(2); e(2) -e(1)];
           T_i  = [1 0; 0 1; -ro/1000*e(2)   ro/1000*e(1)];
           v_Forces(:,i) = T_i*beta_i*[Fn_ij; Fs_ij];
           
       end
        
    end
    
    
    
end


end

