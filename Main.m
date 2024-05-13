close;
clear;
clc;

%% Salt and Pepper Noise removal Using type 2 fuzzy system

input_image=imread('Images/lena.png');
figure(1), imshow(input_image);
title('input image');
[~,~,d] = size(input_image);
if d == 3
    im_gray=rgb2gray(input_image);
    im_gray_1=im2double(im_gray);
else
    im_gray_1=im2double(input_image);
end
figure(2), imshow(im_gray_1);
title('Pillars');

% time=tic;
% % input_image=imread('peppers.bmp');
% input_image=imread('woman.jpg');
% figure(1),imshow(input_image);
% % im_gray=rgb2gray(input_image);
% im_gray=input_image;
% im_gray_1=im2double(im_gray); % Will be easy to calculate PSNR in the end
% figure(2),imshow(im_gray_1);

%% Inialization Of Parameters
nd = 0.90; Salt and pepper nosie density
im_noised=imnoise(im_gray_1,'salt & pepper',nd);
figure(3),imshow(im_noised);
title('Pillars');
[p,q]=size(im_noised);

% % Inialising Denoised image %% Image size is being increased by 8, to take% care of edges.
% im_denoised=0.53*ones(p+16,q+16);
% im_denoised(9:p+8,9:q+8)=im_noised; %
% 
% % M=1;  % M is half window size
% count0=0;
% count1=0;
% M=1;  % M is half window size %% Inialization
% N_init = 8;   % for low noise image N can be low
% im_denoised_pixels = zeros(p+16,q+16);
% im_noised_pixels = zeros(p+16,q+16);

%% new
[p,q]=size(im_noised);
im_denoised=0.73*ones(p+20,q+20);% Padding for edges
im_denoised(10:p+9,10:q+9)=im_noised; %Image at center
count0=0;
count1=0;
M=1;  % M is half window size %% Inialization
N_init = 8;   % for low noise image N can be low

%% Type-2 Fuzzy Identifier
im_denoised_pixels = zeros(p+20,q+20);%denoised image
im_noised_pixels = zeros(p+20,q+20);%mask showing location of SAP


%% time_elapsed_per_itteration = zeros(1,10);
tic
e=1;
for z=1:10 % For loop for Iterations
    im_denoised_pixels = zeros(p+20,q+20);
    im_noised_pixels = zeros(p+20,q+20);
    %     c=dm;  % To terminate the loop (not necessary to use)
    for j=10:q+9 % To scan rows
        for i=10:p+9 % To scan col.
             if (nd <= 0.30)
                M = 1;  % M is half window size %% Inialization 
            elseif (0.30 < nd <= 80)
                M = 2;  % M is half window size %% Inialization 
            else 
                M = 3;  % M is half window size %% Inialization 
            end
            N_init = 6;   % for low noise image N can be low
            %
            S_max = 2; % upper bound of 'M'
            N = N_init; % Inialisation of number of good pixel required to evaluate the pixel value of a currupted pixel

            while (im_denoised(i,j)==0)||(im_denoised(i,j)==1)
               
                PixelVec = PixelVector(im_denoised,i,j,M); % Vector containing pixel around possible currupted pixel
                o = 1;
%% finding of adaptive threshold
                if (nnz(PixelVec==0) & nnz(PixelVec==1));

                    G_1 = PixelVec(PixelVec~=0 & PixelVec~=1);
                    ALD = abs(G_1 - sum(PixelVec(~(PixelVec~=0 & PixelVec~=1)))/length(PixelVec(~(PixelVec~=0 & PixelVec~=1))));
                    if isempty(G_1)
                        M=M+1;
                        continue
                    end
                elseif (nnz(PixelVec==0) & ~nnz(PixelVec==1));
                    G_1 = PixelVec(PixelVec~=0);
                    ALD = abs(G_1 - sum(PixelVec(~(PixelVec~=0)))/length(PixelVec(~(PixelVec~=0))));
                    if isempty(G_1)
                        M=M+1;
                        continue
                    end
                else  (~nnz(PixelVec==0) & nnz(PixelVec==1));
                    G_1 = PixelVec(PixelVec~=1);
                    ALD = abs(G_1 - sum(PixelVec(~(PixelVec~=1)))/length(PixelVec(~(PixelVec~=1))));
                    if isempty(G_1)
                        M=M+1;
                        continue
                    end
                end
                lenght_G_1 = length(G_1);
                MALD = max(ALD);
                [T_min,T_max,T_min_max,PI,H,sigma,average_mu] = T2MF(G_1);
                T_Threshold = T_max;
                if isvector(G_1)
                    ave_PI = PI;
                else
                    ave_PI = sum(PI)/H;
                end
                
                if(MALD <= T_min)
                    L =abs(T_min);

                elseif (T_min<MALD) && (MALD<T_max)

                    L= abs((MALD -T_min)/(T_max-T_min));

                else
                    L = abs(T_max);

                end
%% Selection of good pixels               
                for x=1:lenght_G_1
                    if  ave_PI(x)>=T_Threshold 
                        count0=count0+1;
                        G(o)=G_1(x);
                        o=o+1;
                    elseif  (G_1(x)~=0)&&(G_1(x)~=1)
                        count1=count1+1;
                        G(o)=G_1(x);
                        o=o+1;
                    end
                    if o==N+1
                        break
                    end
                end

                neta=length(find(G));

                if (neta<N) && (M< S_max)
                    M=M+1;

                    continue
                elseif (neta<N) && (M== S_max)
                    N=N-1;
                    if N<1
                        S_max=S_max+1;
                        N=1;
                    end
                    continue
                end
%% denoising of noisy pixels from the good pixels

                % for k=1:(length(G)/2)
                %     Mean(k)=KMM(k,G);
                % end
                % mean_G=((sum(Mean))/length(G));
                % var_G=4.5*abs(G-mean_G);                %% for the Building sigma = 2.5 for rest it is 4.5
                % var_G=mean(var_G);
                % if var_G<=.01;
                %     var_G=.01;
                % end
                mean_G=mean(G);
                [~,~,var_G]=ksdensity(abs(G-mean_G));                
                if var_G<=.01;
                    var_G=.01;
                end
                
                w=gaussmf(G,[var_G,mean_G]);
                W=sum(w);
                weighted_G=w*G';
                im_denoised_pixels(i,j) = (1-L)*im_denoised(i,j) + L*weighted_G/W;
                im_noised_pixels(i,j) = im_denoised(i,j);
                break
            end
        end
    end
    %time_elapsed_per_itteration(z)=toc;
     
    im_denoised=(im_denoised-im_noised_pixels)+im_denoised_pixels;
end
time_elapsed=toc;

%%
im_denoised=im_denoised(10:p+9,10:q+9);
im_denoised_1 = im2uint8(im_denoised);
im_denoised_1=int16(im_denoised_1);
im_gray=int16(im_gray);
PSNR=10*log10((255*255)/((1/((p-10)*(q-10)))*sum(sum((im_denoised_1(10:p-9,10:q-9)-im_gray(10:p-9,10:q-9)).^2))));
% PSNR=10*log10((255*255)/((1/((p-10)*(q-10)))*sum(sum((im_denoised_1(6:p-5,6:q-5)-im_gray(6:p-5,6:q-5)).^2))));
% PSNR=10*log10(1/((1/((p-14)*(q-14)))*sum(sum((im_gray(8:p-7,8:q-7)-im_denoised(8:p-7,8:q-7)).^2))));
% PSNR=10*log10(1/((1/((p)*(q)))*sum(sum((im_gray-im_denoised).^2))));
PSNR_char=num2str(PSNR);
figure(4);
imshow(im_denoised);
% title(PSNR_char);

