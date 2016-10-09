%                            ______________________
%                            
%                            manual 3d image resize
%                            ______________________
%

function out = imresize3d(A,numpix)              % A - input image, 
                                                 % numpix - new number of pixels in each dimension

                            
% A = rand(53,63,46);


% say i want 16x16x16

ratio = size(A)/numpix;

r = ratio;

% NOTE:

% using round(r(1)*i) internally in the loop, this way every 3 out of 10 pixels
% for the 1st dimension (rows) contains an 'extra original pixel'. more
% symmetrical resampling this way

% added in round(r(2)*j) and round(r(3)*k)



% want to take A(1:3,1:3,1:3) and get the mean

A_rs = zeros(numpix,numpix,numpix);  % resized matrix
A_stor = 0;              % matrix to store components and squeeze into vector before summing


endpt = size(A);

for i=1:numpix
    
    
    for j = 1:numpix
                                           % -v- may not be needed if I use the rounding method
                                           
        for k = 1:numpix                   % if statement to manage end points of dimensions
            if (i~=numpix & j~=numpix)          % 2^3 = 8 possibilites to check for
                if (k==numpix)
                    A_stor = A(round((r(1)*(i-1))+1):round((r(1)*i)), round((r(2)*(j-1))+1):round((r(2)*j)), round((r(3)*(k-1))+1):endpt(3));
                    A_rs(i,j,k) = mean(A_stor(:));
                else
                    A_stor = A(round((r(1)*(i-1))+1):round((r(1)*i)), round((r(2)*(j-1))+1):round((r(2)*j)), round((r(3)*(k-1))+1):round((r(3)*k)));
                    A_rs(i,j,k) = mean(A_stor(:));
                end
           
               
            elseif (i~=numpix & j==numpix)
                if (k==numpix)
                    A_stor = A(round((r(1)*(i-1))+1):round((r(1)*i)), round((r(2)*(j-1))+1):endpt(2), round((r(3)*(k-1))+1):endpt(3));
                    A_rs(i,j,k) = mean(A_stor(:));
                else
                    A_stor = A(round((r(1)*(i-1))+1):round((r(1)*i)), round((r(2)*(j-1))+1):endpt(2), round((r(3)*(k-1))+1):round((r(3)*k)));
                    A_rs(i,j,k) = mean(A_stor(:));
                end
            
            
            elseif (i==numpix & j~=numpix)
                if (k==numpix)
                    A_stor = A(round((r(1)*(i-1))+1):endpt(1), round((r(2)*(j-1))+1):round((r(2)*j)), round((r(3)*(k-1))+1):endpt(3));
                    A_rs(i,j,k) = mean(A_stor(:));
                else
                    A_stor = A(round((r(1)*(i-1))+1):endpt(1), round((r(2)*(j-1))+1):round((r(2)*j)), round((r(3)*(k-1))+1):round((r(3)*k)));
                    A_rs(i,j,k) = mean(A_stor(:));
                end
            
            
            elseif (i==numpix & j==numpix)
                if (k==numpix)
                    A_stor = A(round((r(1)*(i-1))+1):endpt(1), round((r(2)*(j-1))+1):endpt(2), round((r(3)*(k-1))+1):endpt(3));
                    A_rs(i,j,k) = mean(A_stor(:));
                else
                    A_stor = A(round((r(1)*(i-1))+1):endpt(1), round((r(2)*(j-1))+1):endpt(2), round((r(3)*(k-1))+1):round((r(3)*k)));
                    A_rs(i,j,k) = mean(A_stor(:));
                end
            end
          
            
        end
        
    end
    
end


out = A_rs;


% figure; 
% subplot(1,2,1); imagesc(A(:,:,1), [0 1]); colormap 'gray';
% subplot(1,2,2); imagesc(A_rs(:,:,1), [0 1]); colormap 'gray';
% 


% sum(A(:))                % to check if count of A and A_rs is same change 'mean' to 'sum' in for loops
% sum(A_rs(:))

