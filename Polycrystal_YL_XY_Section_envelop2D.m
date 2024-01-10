clear;      
close;

%% Creating Bishop Hill stress state matrix from the text file named: BHfile.txt

B = fopen('BHfile.txt');
BH = textscan(B, ' %f %f %f %f %f %f');
fclose(B);

% %% ppt details
% f = 0.008;
% sigma_bar = 10000e6;
% tau = 88e6;

%% Reading the orientation file

prompt = 'The euler angle file path and name with .txt extension \n';
g_vectorfile = input(prompt);                        
g = fopen(g_vectorfile);                            
g_matrix = textscan(g, '%f %f %f');               
l_g =  length(g_matrix{1,1});

%% Reading the strain file

S = fopen('strains.txt');
strain = textscan(S, ' %f %f %f ');
l_s =  length(strain{1,1});
fclose(S);
M = zeros(1,l_s);
ro = zeros(1,l_s);


     for u=1:1:l_s
         e_ext=[strain{1,1}(u),0,0;0,strain{1,2}(u),0;0,0,strain{1,3}(u)];
         Wmax= zeros(l_g,1);
         
         for c=1:1:l_g                     
            A = DC_matrix_function(g_matrix{1,1}(c),g_matrix{1,2}(c),g_matrix{1,3}(c));
           [e]= transform_e_function(e_ext,A);
            W= zeros(1,56);
            BH_state = zeros(56,6);
            
                    for m=1:1:56 
                        W(m)= -(BH{1,2}(m)*e(1,1))+ BH{1,1}(m)*e(2,2)+ BH{1,4}(m)*(e(2,3)+e(3,2))+BH{1,5}(m)*(e(1,3)+e(3,1))+BH{1,6}(m)*(e(1,2)+e(2,1));
                        BH_state(m,:) = [BH{1,1}(m),BH{1,2}(m),BH{1,3}(m),BH{1,4}(m),BH{1,5}(m),BH{1,6}(m)]; % [A,B,C,F,G,H]
                    end                    
                    Wmax(c)= max(abs(W));        
        end
        M(u) = mean(Wmax)/e_ext(1,1);
        ro(u) = -strain{1,2}(u)/strain{1,1}(u); 
     end 

%% Matix to store YL data points
th = 0:1:360;
radius = 2*ones(size(th));

for u=1:1:l_s    
        %aX+bY = M, a =1; b = -rho       
        b= -ro(u);
        m = M(u);
        r = abs(m./(cosd(th) + b.*sind(th)));
        radius = bsxfun(@min,r,radius);        
end   

%% Plotting the YL
[X,Y] = pol2cart(th*(pi/180),radius);
plot(X,Y)
pbaspect([1 1 1])
xlim([-2 2]);
ylim([-2 2]);
xlabel('SigmaXX');
ylabel('SigmaYY');
