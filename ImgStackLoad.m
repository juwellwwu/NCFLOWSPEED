function [ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel)

fprintf('Loading input images...\n');

if InputImg_TifvsMAT_Query==1    
    
        ImgSlices=dir(fullfile(BatchImgInputFolder,'*.tif'));
        NumImgSlices=length(ImgSlices); 
        
        % Get ImgSize from 1st file 
        ImgNameString=fullfile(BatchImgInputFolder,ImgSlices(1).name);
        [Img,~,~]=imread(ImgNameString);
        ImgStack=zeros(size(Img,1),size(Img,2),NumImgSlices,'single');
        
        for k = 1:NumImgSlices   
            ImgNameString=fullfile(BatchImgInputFolder,ImgSlices(k).name);
            [Img,~,~]=imread(ImgNameString);  % Img: class uint8; 3 to 4 layers (R-G-B-alpha)
            
            % if RGB image, split and choose channel
            if size(Img,3)>=3
                if RGBChannel==1
                    Img=Img(:,:,1);
                elseif RGBChannel==2
                    Img=Img(:,:,2);
                elseif RGBChannel==3
                    Img=Img(:,:,3);
                end
            end
            
% % %             if k==1       
% % %                     ImgStack=zeros(size(Img,1),size(Img,2),NumImgSlices,'single');
% % %             end           
            ImgStack(:,:,k)=Img(:,:,1);      
        end
        
        % 17Feb27: Single precision float calculation
        ImgStack=single(mat2gray(ImgStack,[0,255]));
         
elseif InputImg_TifvsMAT_Query==2
    
        ImgStack=single(importdata(ImgStackData_FilenameString)./255);
                        
end

Img_Height=size(ImgStack,1);
Img_Width=size(ImgStack,2);
NumImgSlices=size(ImgStack,3);

end