%% This code calculates Yield Locus for a single crystal. It requires only 
% single orientation as the input text file. This code calculates stress
% states active in the orientations and from it it takes the sigma xx and
% yy to plot the yiled locus.
% This script uses BH method to calculate the YL for a single orientation.
% This works for fcc materials

clear;      
close;

%% Creating Bishop Hill stress state matrix from the text file named: BHfile.txt

B = fopen('BHfile.txt');
BH = textscan(B, ' %f %f %f %f %f %f');
fclose(B);

%% Reading the orientation file

prompt = ['The euler angle file name with .txt extension \n'...
            'You can chose from the predefined orientaions\n' ...
            ' cube, Q, R, s,goss, brass etc, or make your own\n'...
            ' similar to these\n'];
        
g_vectorfile = input(prompt);                        
g = fopen(g_vectorfile);                            
g_matrix = textscan(g, '%f %f %f');               
l_g =  length(g_matrix{1,1});

%% Reading the strain file

S = fopen('strains.txt');
strain = textscan(S, ' %f %f %f ');
l_s =  length(strain{1,1});
fclose(S);
YieldLocusData = struct();
CrystalData =struct();

    for c=1:1:l_g
%         
        for u=1:1:l_s
            
            e_ext=[strain{1,1}(u),0,0;0,strain{1,2}(u),0;0,0,strain{1,3}(u)];
            A = DC_matrix_function(g_matrix{1,1}(c),g_matrix{1,2}(c),g_matrix{1,3}(c));
            [e]= transform_e_function(e_ext,A);
            W= zeros(56,1);
            BH_state = zeros(56,6);
            
                    for m=1:1:56 

                        W(m)= -(BH{1,2}(m)*e(1,1))+ BH{1,1}(m)*e(2,2)+ BH{1,4}(m)*(e(2,3)+e(3,2))+BH{1,5}(m)*(e(1,3)+e(3,1))+BH{1,6}(m)*(e(1,2)+e(2,1));
                        BH_state(m,:) = [BH{1,1}(m),BH{1,2}(m),BH{1,3}(m),BH{1,4}(m),BH{1,5}(m),BH{1,6}(m)]; % [A,B,C,F,G,H]
                    end
                    
                    Wmax= max(abs(W));
                    p = 1;                   
                     % Selecting the appropriate Bishop Hill stress state, i.e. which maximises the virtual work                    
                    for n=1:1:56                       
                       
                        if (abs(W(n))== Wmax)

                                                   
                            % Now calculating the appropriate stresses in xtal ref
                            % frame from the bishop hill constansts A,B,C,
                            % F,G,H
                            
                            s23=BH_state(n,4);
                            s31=BH_state(n,5);
                            s12=BH_state(n,6);

                            s11= -(2*(BH_state(n,6)*A(1,3)*A(2,3)+BH_state(n,5)*A(1,3)*A(3,3)+...
                                BH_state(n,4)*A(2,3)*A(3,3)) + BH_state(n,1)*A(2,3)^2 + ...
                                BH_state(n,2)*(A(2,3)^2 + A(3,3)^2))/(A(1,3)^2 + A(2,3)^2 + A(3,3)^2);
                            s22= s11 + BH_state(n,1)+BH_state(n,2);
                            s33= s11 + BH_state(n,2);
                            s13=s31;
                            s32=s23;
                            s21=s12;

                            s_xtal= [s11, s12, s13; s21, s22, s23; s31, s32, s33]; 
                            
                            % Now convert the stresses from xtal coordinates to external reference frame 

                           s = [0,0,0;0,0,0;0,0,0];

                                           % Now calculating transformed stress tensor in external ref frame

                                           for i=1:1:3
                                               for j=1:1:3
                                                   for k=1:1:3
                                                       for l=1:1:3

                                                           s(i,j)= s(i,j)+ A(k,i)*A(l,j)*s_xtal(k,l);

                                                       end
                                                   end
                                               end
                                           end
                                stress(p,:) = [s(1,1),s(2,2),s(3,3),s(2,3),s(1,3),s(1,2)];
                                p= p+1;                                                           

                        end
                    end
        ro = strain{1,2}(u)/strain{1,1}(u); 
        CrystalData(u).Stress_set = stress;
        CrystalData(u).Rho = ro;
        YieldLocusData(c).Ori= [g_matrix{1,1}(c),g_matrix{1,2}(c),g_matrix{1,3}(c)];
        YieldLocusData(c).Data = CrystalData;
 
        end
        
    end 
     
save('YieldLocus_Data')

sigma_xx = zeros(1,55);
sigma_yy = zeros(1,55);
slope = zeros(1,55);
Data = YieldLocusData;
p = 1;
rho_l = length(Data(1,1).Data);
for i = 1:1:55
for j = 1:1:2
sigma_xx(p) = Data.Data(1,i).Stress_set(j,1);% j=1
sigma_yy(p) = Data.Data(1,i).Stress_set(j,2);
slope(p)= Data.Data(1,i).Rho;
p = p+1;
end
end

%%
for k = 1:1:(p-1)
    x1 = sigma_xx(k);
    y1 = sigma_yy(k);
    X = (x1-3):0.1:(x1+3);
    m = -1/slope(k);
    if (slope(k) == 0) 
        Y = -2:0.2:2;
        X = x1*ones(numel(Y));
    else
    Y = m*(X - x1)+y1;
    
    end       
    plot(X,Y,'-')
    hold on
end

h = plot(sigma_xx, sigma_yy, '+');
xlim([-2 2]);
ylim([-2 2]);
% set(h,'XLim',[-2, 2],'YLim',[-2, 2]);
hold off


