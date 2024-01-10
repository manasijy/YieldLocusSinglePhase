% clear;      
% close;

%% Creating Bishop Hill stress state matrix from the text file named: BHfile.txt

B = fopen('BHfile.txt');
BH = textscan(B, ' %f %f %f %f %f %f');
fclose(B);

%% Reading the orientation file

prompt = 'The euler angle file name with .txt extension \n';
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
         for v = 1:1:11
             gamma12 = -1 + 0.2*(v-1);
         e_ext=[strain{1,1}(u),gamma12,0;gamma12,strain{1,2}(u),0;0,0,strain{1,3}(u)];
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
        M(u,v) = mean(Wmax)/e_ext(1,1);
        ro(u,v) = -strain{1,2}(u)/strain{1,1}(u);
        g(u,v) = gamma12/strain{1,1}(u);
         end
     end 

%% Matix to store YL data points
azimuth = 0:5:360;
elevation = -90:5:90;
[az, el] = meshgrid(azimuth, elevation);
radius = 2*ones(size(az));

for u=1:1:l_s
    for v=1:1:11
        %aX+bY+cZ = M, a =1; b = -rho,c = 2gamma        
        b= -ro(u,v);
        c= 2*g(u,v);
        m = M(u,v);
        r = abs(m./(cosd(az).*cosd(el)+ b.*cosd(el).*sind(az)+ c.*sind(el)));
        radius = bsxfun(@min,r,radius);        
    end
end

%% Plotting the YL
[X,Y,Z] = sph2cart(az*(pi/180),el*(pi/180),radius);
surf(X,Y,Z)
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
zlim([-1.5 1.5]);
xlabel('SigmaXX');
ylabel('SigmaYY');
zlabel('TauXY');
