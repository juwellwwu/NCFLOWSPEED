clc; fclose('all'); clearvars;
close all hidden; 
 
%% =====DESCRIPTION=====

% Generates vessel blood flow velocity map

% == Usage: 
% User specifies variables in "USER INPUT" section.

% == Output folders:
% (Each folder includes an 8-bit image stack)
% "HL_Mask": Hough line position in INPUT time series, named by index (-XXXX) 
% "xDistyTime_LinFitSlopeQual": Hough line slope estimate by Radon
% transform and quality of estimate, overlaid on Distance-Time image of each vessel segment

% == Output files:
% "xDyT_SlopeSpeed.mat": vessel segment flow velocities by Hough line slope measurements; quality of measurements
% "Abs FlowSpeed Map1 r300.tif": blood flow velocity map
% "Speed_Map_abs_flipud.mat": blood flow velocity map numerical data
% "Speed_HLIdx_Map_flipud.mat": map w/ index of Hough Line used to create blood flow velocity map data


%%  =====DO NOT REMOVE=====

% Supplementary software code for Jung et al. "Intravital fluorescence microscopy with negative contrast"
% Author: Juwell W. Wu 
% Wellman Center for Photomedicine, Massachusetts General Hospital, Harvard Medical School, Boston, MA 02114, USA 
% Email address: jwwu@@mgh.harvard.edu  
% Last revision: May-2021


%% USER INPUT

addpath(genpath('OUTPUT/')) % Add OUTPUT folder and subfolders to path

% === INPUT: Time Series of Vessel Blood Flow
InputImg_TifvsMAT_Query=1; % "1" for TIF; "2" for MAT

% TIF image directory (do NOT include / at end), if InputImg_TifvsMAT_Query=1
BatchImgInputFolder='INPUT/ImgStack'; 

% MAT file name, if InputImg_TifvsMAT_Query=2
ImgStackData_FilenameString='';

% === INPUT: Output folder name (do NOT include / at end)
BatchImgOutputFolder='OUTPUT/NCFlowSpeed_Out'; 

% === INPUT: Pixel length in XY (um), in Z (frame rate in seconds)
XY_PxLength=0.77;  
Z_PxLength=1/120; 
  
% === INPUT: Segment Length (um)
Min_HL_Length=15;   
Max_HL_Length=60;
 
% === INPUT: Vessel diameter (um)
SmallestVesselSize=4;   
LargestVesselSize=16;  
 
 
%% Optimized inputs: modify with care

% === INPUT: Color of ImgStack
% Input image RGB Channel Number (1=R;2=G;3=B)
% Applicable only for RGB TIF input
RGBChannel=2;


%% Load Image Stack
 
% Load image stack
[ImgStack,Img_Height,Img_Width,NumImgSlices]=ImgStackLoad(InputImg_TifvsMAT_Query,BatchImgInputFolder,ImgStackData_FilenameString,RGBChannel);
 
 
%% METHOD1: Sum Stacks to create ZMask
 
fprintf('Summing Slices...\n');
 
ImgStack_Input=ImgStack;
 
Img_StackSum=sum(ImgStack_Input,3);
 
% Rescale matrix such that lowest val=0; highest val=1
Img_StackSum=(Img_StackSum-min(Img_StackSum(:)))*1/(max(Img_StackSum(:))-min(Img_StackSum(:)));
 
clearvars *_Input;
 
  
%% Noise Removal, StackSum
 
fprintf('Noise Removal...\n');
 
Img_Input=Img_StackSum;
 
strel_Radius=round(SmallestVesselSize/XY_PxLength*1.0);  
SE=strel('disk',strel_Radius);
MedianOrd=round(numel(find(SE.Neighborhood))*0.5);
 
Img_NR=Img_Input;
Img_NR_Med=ordfilt2(Img_NR,MedianOrd,SE.Neighborhood,'symmetric'); 
Img_NR_Idx=find(Img_NR-Img_NR_Med>0.001);  
Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
 
Img_NR_Med=ordfilt2(Img_NR,MedianOrd,SE.Neighborhood,'symmetric');
Img_NR_Idx=find(Img_NR_Med-Img_NR>0.001);   
Img_NR(Img_NR_Idx)=Img_NR_Med(Img_NR_Idx);
 
clearvars Img_StackZMax;
clearvars Img_NR* -except Img_NR;
clearvars MedianOrd strel_Radius SE;
clearvars *_Input;
 
 
%% Image Opening: Erosion followed by Dilation - Remove small objects
 
fprintf('Image opening (erosion, followed by dilation) ...\n');
 
Img_Input=Img_NR;
 
strel_Radius=round(SmallestVesselSize/XY_PxLength-1);  
SE=strel('disk',strel_Radius);
 
Img_Min=imerode(Img_Input,SE);
 
Img_Input=Img_Min;
Img_Max=imdilate(Img_Input,SE);
 
Img_Open=single((Img_Max-min(Img_Max(:)))*1/(max(Img_Max(:))-min(Img_Max(:))));
 
clearvars Img_NR;
clearvars strel_Radius SE;
clearvars Img_Min Img_Max;
clearvars *_Input;
 
 
%% CLAHE
 
fprintf('CLAHE ...\n');
 
Img_Input=Img_Open;
 
CLAHE_RowSize=round(SmallestVesselSize/XY_PxLength*0.5);    
CLAHE_ColSize=CLAHE_RowSize;
Img_CLAHE=adapthisteq(Img_Input,'clipLimit',0.025,'NumTiles',[CLAHE_RowSize CLAHE_ColSize],'Distribution','uniform');

clearvars Img_Open;
clearvars CLAHE_RowSize CLAHE_ColSize;
clearvars *_Input;
 
 
%% Threshold
 
fprintf('Thresholding...\n');
 
Img_Input=Img_CLAHE;
 
Img_Thresh=Img_Input;
 
Img_ThreshLevel=multithresh(Img_Thresh,3); 
 
Img_Thresh=imbinarize(Img_Thresh,Img_ThreshLevel(1,2));
 
clearvars Img_CLAHE;
clearvars Img_ThreshLevel;
clearvars *_Input;
 
 
%% Blob Property Analysis
 
fprintf('Exclude Dim, Small areas from Flow Speed Analysis...\n');
 
Img_Input=Img_Thresh;
Img_StackSum_Input=Img_StackSum;
 
[Img_Blob_Lb, Img_Blob_NumRegions]=bwlabel(Img_Input);
Blob_regionprops=regionprops(Img_Blob_Lb,'PixelIdxList','MajorAxisLength');
 
Blob_MajorAxisLength_List=cat(1,Blob_regionprops.MajorAxisLength);
        
Blob_MajorAxisLength_Flag=bsxfun(@lt,Blob_MajorAxisLength_List,Min_HL_Length/XY_PxLength);
Blob_MajorAxisLength_FlagID=find(Blob_MajorAxisLength_Flag);
 
Blob_MaxIntensity_List=zeros(size(Blob_regionprops,1),1);
for k=1:size(Blob_regionprops,1)
    Mask=zeros(size(Img_Blob_Lb),'logical');
    Mask(Img_Blob_Lb==k)=1;
    StackSum_Masked=Img_StackSum_Input.*Mask;
    Blob_MaxIntensity_List(k,1)=max(StackSum_Masked(:));
end
Blob_MaxIntensity_Flag=bsxfun(@lt,Blob_MaxIntensity_List,0.05);
Blob_MaxIntensity_FlagID=find(Blob_MaxIntensity_Flag);
 
Blob_Prop_FlagID=union(Blob_MajorAxisLength_FlagID,Blob_MaxIntensity_FlagID);
 
Img_PreSkel=Img_Input;
 
ZeroPxLinIdxList=unique(cat(1,Blob_regionprops(Blob_Prop_FlagID).PixelIdxList));
 
Img_PreSkel(ZeroPxLinIdxList)=0;


clearvars Img_Thresh;
clearvars Img_StackSum;
clearvars Mask StackSum_Masked;
clearvars Img_Blob_Lb Img_Blob_NumRegions;
clearvars Blob_regionprops;
clearvars Blob_MajorAxisLength_Flag Blob_MajorAxisLength_FlagID;
clearvars Blob_MaxIntensity_Flag Blob_MaxIntensity_FlagID;
clearvars Blob_MajorAxisLength_List Blob_MaxIntensity_List;
clearvars Blob_Prop_FlagID;
clearvars ZeroPxLinIdxList;
clearvars *_Input;
 
 
%% Skeletonization of Z-Mask of Threshold
 
fprintf('Skeletonizing Z-Mask of Threshold...\n');
 
Img_Input=Img_PreSkel;
 
Skeleton_Length=round(LargestVesselSize/XY_PxLength*2.0);   
Img_Skeleton=logical(bwmorph(skeleton(Img_Input)>Skeleton_Length,'skel',Inf)); 

clearvars Skeleton_Length;
clearvars *_Input;
 
 
%% Widen Skeleton lines
 
fprintf('Maximum Filtering to Widen Skeleton Lines...\n');
 
Img_Input=Img_Skeleton;
 
strel_Size=round(SmallestVesselSize/XY_PxLength-1);  
SE=strel('square',strel_Size);
 
Img_WideSkel=imdilate(Img_Input,SE);
 
clearvars Img_Skeleton;
clearvars strel_Size SE;
clearvars *_Input;
 
 
%% Prepare for saving
 
SaveFilePath=strcat(BatchImgOutputFolder,'/'); 
mkdir(SaveFilePath);
 
 
%% Hough Transform Z-Union Threshold Skeleton
 
fprintf('Hough Transform...\n');
 
Img_WideSkel_Input=Img_WideSkel;
 
[Hough_Matrix,Hough_THETA,Hough_RHO]=hough(Img_WideSkel_Input,'RhoResolution',0.5,'Theta',-90:0.5:89.5); 

HP_NHoodSize=floor(SmallestVesselSize/XY_PxLength*1.0); 
if rem(HP_NHoodSize,2)==0
    HP_NHoodSize=HP_NHoodSize-1;
end
Hough_Peaks=houghpeaks(Hough_Matrix,1E9,'threshold',1,'NHoodSize',[HP_NHoodSize HP_NHoodSize]);
 
Min_HL_Length_Px=round(Min_HL_Length/XY_PxLength);
Hough_Lines_PreFlagRmv1=houghlines(Img_WideSkel_Input,Hough_THETA,Hough_RHO,Hough_Peaks,'MinLength',Min_HL_Length_Px);
 
clearvars Hough_Matrix Hough_THETA Hough_RHO;
clearvars HP_NHoodSize;
clearvars Hough_Peaks;
clearvars Min_HL_Length_Px;
clearvars *_Input;
 
 
%% Create Hough Line Mask, Keep Short Lines with connected components only

fprintf('Evaluating Hough Lines for length...\n');
 
HL_Input=Hough_Lines_PreFlagRmv1;
 
xy_Endpts=[cat(1,HL_Input.point1) cat(1,HL_Input.point2)];
 
Flag_HL_Long=zeros(length(HL_Input),1);
 
HL_PreFlagRmv1_Length=sqrt((xy_Endpts(:,3)-xy_Endpts(:,1)).^2+(xy_Endpts(:,4)-xy_Endpts(:,2)).^2);
 
Flag_HL_Long(HL_PreFlagRmv1_Length>(Max_HL_Length/XY_PxLength))=1;
 
Flag_HL_Long_ldx=find(Flag_HL_Long);

clearvars xy_Endpts;
clearvars HL_PreFlagRmv1_Length;
clearvars Flag_HL_Long;
clearvars *_Input;
 
 
%% Remove Flagged Hough Lines
 
fprintf('Removing Hough Lines Flagged for length ... \n');
 
HL_PreFlagRmv_Input=Hough_Lines_PreFlagRmv1;
Flag_HL_Long_ldx_Input=Flag_HL_Long_ldx;

Hough_Lines_PreFlagRmv2=HL_PreFlagRmv_Input;
 
Hough_Lines_PreFlagRmv2(Flag_HL_Long_ldx_Input)=[];
 
fprintf(strcat(num2str(length(Flag_HL_Long_ldx_Input)),'/',num2str(length(HL_PreFlagRmv_Input)),' Hough Lines removed for length. \n\n'));

clearvars Hough_Lines_PreFlagRmv1;
clearvars Flag_HL_Long_ldx;
clearvars *_Input;
 
 
%% Create Hough Line Mask, Keep Hough Lines with connected components only
 
fprintf('Evaluating Hough Lines for connected components; creating mask ... \n');
 
HL_Input=Hough_Lines_PreFlagRmv2;
 
% === Create endpoint coordinate matrix
% # Rows = # Hough Lines
% Col1 = x-coordinate, Endpoint1
% Col2 = y-coordinate, Endpoint1
% Col3 = x-coordinate, Endpoint2
% Col4 = y-coordinate, Endpoint2
xy_Endpts=[cat(1,HL_Input.point1) cat(1,HL_Input.point2)];
 
xinterp_Query=~(abs(xy_Endpts(:,1)-xy_Endpts(:,3))>=abs(xy_Endpts(:,2)-xy_Endpts(:,4)));
 
Flag_HL_Unconn=zeros(length(HL_Input),1);
 
HL_xDistyTime_PreFlagRmv_Cell=cell(length(HL_Input),1);
 
HL_LinIdx_PreFlagRmv_Cell=cell(length(HL_Input),1);
 
HL_Mask_PreFlagRmv_Stack=zeros(Img_Height,Img_Width,length(HL_Input));
 
for k = 1:length(HL_Input)

    if xinterp_Query(k)==0 
        A1_Endpts=[xy_Endpts(k,1) xy_Endpts(k,3)];
        A2_Endpts=[xy_Endpts(k,2) xy_Endpts(k,4)];
    else        
        A1_Endpts=[xy_Endpts(k,2) xy_Endpts(k,4)];  
        A2_Endpts=[xy_Endpts(k,1) xy_Endpts(k,3)];  
    end
        
    if A1_Endpts(2)>A1_Endpts(1)
        A1_Vector=(A1_Endpts(1):A1_Endpts(2))';  
    else
        A1_Vector=(A1_Endpts(2):A1_Endpts(1))'; 
    end
    A2_Vector=interp1(A1_Endpts,A2_Endpts,A1_Vector);  
    A2_Vector=round(A2_Vector);
    
    if xinterp_Query(k)==0 
        HL_LinIdx=sub2ind(size(ImgStack(:,:,1)),A2_Vector, A1_Vector);            
    else 
        HL_LinIdx=sub2ind(size(ImgStack(:,:,1)),A1_Vector, A2_Vector);  
    end
       
    HL_Mask=zeros(size(ImgStack(:,:,1)),'logical');
    HL_Mask(HL_LinIdx)=1;
    
    HL_ConnTest=HL_Mask.*Img_WideSkel; 
    HL_ConnTest(~HL_Mask)=1;  
    if length(find(~HL_ConnTest))/length(find(HL_Mask))>0.1    
        Flag_HL_Unconn(k,1)=1;
    end
          
    HL_IntensityList=ImgStack(find(repmat(HL_Mask,1,1,size(ImgStack,3))));
    HL_xDistyTime=rot90(reshape(HL_IntensityList,size(A1_Vector,1),NumImgSlices),3);
 
    HL_xDistyTime_PreFlagRmv_Cell{k,1}=single(HL_xDistyTime);
    HL_LinIdx_PreFlagRmv_Cell{k,1}=HL_LinIdx;
    HL_Mask_PreFlagRmv_Stack(:,:,k)=HL_Mask; 
            
end
 
Flag_HL_Unconn_Idx=find(Flag_HL_Unconn);

clearvars HL_LinIdx;
clearvars HL_Mask;
clearvars HL_ConnTest;
clearvars HL_IntensityList HL_xDistyTime;
clearvars A1_* A2_*;
clearvars xy_Endpts;
clearvars xinterp_Query;
clearvars Flag_HL_Unconn;
clearvars *_Input;
 
 
%% Remove Flagged Hough Lines
 
fprintf('Removing Flagged Hough Lines (for length, unconnected components)...\n');
 
HL_Input=Hough_Lines_PreFlagRmv2;
xDyT_Cell_Input=HL_xDistyTime_PreFlagRmv_Cell;
Linldx_CeIl_Input=HL_LinIdx_PreFlagRmv_Cell;
Mask_Stack_Input=HL_Mask_PreFlagRmv_Stack;
 
Hough_Lines=HL_Input;
HL_xDistyTime_Cell=xDyT_Cell_Input;
HL_LinIdx_Cell=Linldx_CeIl_Input;
HL_Mask_Stack=Mask_Stack_Input; 
 
Hough_Lines(Flag_HL_Unconn_Idx)=[];
HL_xDistyTime_Cell(Flag_HL_Unconn_Idx)=[];
HL_LinIdx_Cell(Flag_HL_Unconn_Idx)=[];
HL_Mask_Stack(:,:,Flag_HL_Unconn_Idx)=[];
 
fprintf(strcat(num2str(length(Flag_HL_Unconn_Idx)),'/', num2str(length(HL_Input)),' Hough Lines removed for unconnected components.\n'));

clearvars Hough_Lines_PreFlagRmv2 HL_xDistyTime_PreFlagRmv_Cell HL_LinIdx_PreFlagRmv_Cell HL_Mask_PreFlagRmv_Stack;
clearvars Flag_HL_Unconn_Idx;
clearvars HL_LinIdx_Cell;
clearvars *_Input;
 
 
%% Display Hough Lines (Short, Connected Components Only)
 
fprintf('Displaying Hough Lines (Short, Connected Components Only)...\n');
 
HL_Input=Hough_Lines;
Img_PreSkel_Input=Img_PreSkel;
Img_WideSkel_Input=Img_WideSkel;
 
figHandle05=figure; 
imshow(Img_PreSkel_Input); 
hold on;
 
for k = 1:length(HL_Input)

   xy=[Hough_Lines(k).point1; Hough_Lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',0.5,'Color','blue');
 
   plot(xy(1,1),xy(1,2),'x','LineWidth',0.5,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',0.5,'Color','red');
   
end
 
% print(figHandle05,'-dtiffn','-r0',strcat(SaveFilePath,'HoughLines_PreSkel.tif'));
 
figHandle06=figure; 
imshow(Img_WideSkel_Input); 
hold on;
 
for k = 1:length(HL_Input)
   xy=[Hough_Lines(k).point1; Hough_Lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',0.5,'Color','blue');
 
   plot(xy(1,1),xy(1,2),'x','LineWidth',0.5,'Color','yellow');    
   plot(xy(2,1),xy(2,2),'x','LineWidth',0.5,'Color','red'); 
end
 
% print(figHandle06,'-dtiffn','-r0',strcat(SaveFilePath,'HoughLines_Skeleton.tif'));
 
clearvars Hough_Lines;
clearvars Img_WideSkel;
clearvars xy;
clearvars *_Input;
 
 
%% Save HL_Mask images
 
fprintf('Saving Hough Line Mask images...\n');
 
HL_Mask_Stack_Input=HL_Mask_Stack;
 
close all;
 
SaveHLMaskFilePath=strcat(SaveFilePath,'HL_Mask/');
mkdir(SaveHLMaskFilePath);
 
for k=1:size(HL_Mask_Stack_Input,3)
   imwrite(double(HL_Mask_Stack_Input(:,:,k)),strcat(SaveHLMaskFilePath,'HLMask -',num2str(k,'%04.0f'),'.tif'));
end
 
clearvars *_Input;
 
 
%% Clear noise in xDistyTime images
 
fprintf('xDistyTime Image noise cleanup...\n');

% SavexDyTBCAdj0MeanFilePath=strcat(SaveFilePath,'xDyT_BCAdj0Mean/');
% mkdir(SavexDyTBCAdj0MeanFilePath);
 
xDyT_Cell_Input=HL_xDistyTime_Cell;
 
HL_xDistyTime_BCAdj0Mean_Cell=cell(size(xDyT_Cell_Input,1),1);
 
for k=1:size(xDyT_Cell_Input,1)         
    
    xDyT=xDyT_Cell_Input{k,1};
    
    xDyT_ColMean=mean(xDyT,1,'omitnan');
    xDyT_ColStdev=std(xDyT,0,1,'omitnan');
    xDyT_ColStdev_div_Mean=xDyT_ColStdev./xDyT_ColMean;
    xDyT_ColStdev_ColFlagIdx=find(xDyT_ColStdev_div_Mean<max(xDyT_ColStdev_div_Mean(:))*0.8);  
    
    xDyT_RowMean=mean(xDyT,2,'omitnan');
    xDyT_RowStdev=std(xDyT,0,2,'omitnan');
    xDyT_RowStdev_div_Mean=xDyT_RowStdev./xDyT_RowMean;
    xDyT_RowStdev_RowFlagIdx=find(xDyT_RowStdev_div_Mean<max(xDyT_RowStdev_div_Mean(:))*0.25); 
    
    xDyT_Temp=xDyT;
    xDyT_Temp(:,xDyT_ColStdev_ColFlagIdx)=NaN;
    xDyT_Temp(xDyT_RowStdev_RowFlagIdx,:)=NaN;
    xDyT_Temp=reshape(xDyT_Temp,numel(xDyT_Temp),1);
    xDyT_Temp(isnan(xDyT_Temp))=[];
    NumHistBins=round(numel(xDyT_Temp)/10);
    HistCtEdge=min(xDyT_Temp(:)):1/NumHistBins:max(xDyT_Temp(:));
    HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
    [HistCt,~]=histcounts(xDyT_Temp,HistCtEdge);
    [~,BkgdHistBinIdx]=min(abs(cumsum(HistCt)-(0.1E-2)*(numel(xDyT_Temp)))); 
    xDyT_Bkgd=HistCt_BinLoc(BkgdHistBinIdx);    
    xDyT_Bkgd_FlagIdx=find(xDyT<xDyT_Bkgd);

    xDyT_BCAdj=xDyT; 
    xDyT_BCAdj=(xDyT_BCAdj-xDyT_Bkgd)*1.0/(max(xDyT_Temp(:))-xDyT_Bkgd);
    xDyT_BCAdj(xDyT_Bkgd_FlagIdx)=NaN;
    xDyT_BCAdj(xDyT_BCAdj<0)=NaN;
    
    xDyT_Temp2=xDyT_BCAdj;
    xDyT_Temp2_ColMean_ColFlagIdx=[];
    xDyT_Temp2_ColMean=mean(xDyT_Temp2,1,'omitnan');
    for ii=1:3        
        xDyT_Temp2_ColMean_Mean=mean(xDyT_Temp2_ColMean,2,'omitnan');
        xDyT_Temp2_ColMean_Stdev=std(xDyT_Temp2_ColMean,0,2,'omitnan');
        xDyT_Temp2_ColMean_ColFlagIdx=...
            intersect(find(xDyT_Temp2_ColMean>=xDyT_Temp2_ColMean_Mean+2*xDyT_Temp2_ColMean_Stdev),...
            find(xDyT_Temp2_ColMean<=xDyT_Temp2_ColMean_Mean-2*xDyT_Temp2_ColMean_Stdev));
        xDyT_Temp2_ColMean(1,xDyT_Temp2_ColMean_ColFlagIdx)=NaN;
    end
    
    xDyT_Temp2=xDyT_BCAdj;
    xDyT_Temp2(:,xDyT_Temp2_ColMean_ColFlagIdx)=NaN;
    xDyT_Temp2_RowMean_RowFlagIdx=[];
    xDyT_Temp2_RowMean=mean(xDyT_Temp2,2,'omitnan');
    for ii=1:3        
        xDyT_Temp2_RowMean_Mean=mean(xDyT_Temp2_RowMean,1,'omitnan');
        xDyT_Temp2_RowMean_Stdev=std(xDyT_Temp2_RowMean,0,1,'omitnan');
        xDyT_Temp2_RowMean_RowFlagIdx=...
            intersect(find(xDyT_Temp2_RowMean>=xDyT_Temp2_RowMean_Mean+2*xDyT_Temp2_RowMean_Stdev),...
            find(xDyT_Temp2_RowMean<=xDyT_Temp2_RowMean_Mean-2*xDyT_Temp2_RowMean_Stdev));
        xDyT_Temp2_RowMean(xDyT_Temp2_RowMean_RowFlagIdx,1)=NaN;
    end
   
    ColFlagIdx=unique(union(xDyT_ColStdev_ColFlagIdx,xDyT_Temp2_ColMean_ColFlagIdx))
    RowFlagIdx=unique(union(xDyT_RowStdev_RowFlagIdx,xDyT_Temp2_RowMean_RowFlagIdx))
    
    xDyT_BCAdj2=xDyT_BCAdj;
    xDyT_BCAdj2=xDyT_BCAdj2./repmat(xDyT_Temp2_ColMean,size(xDyT,1),1);
    xDyT_BCAdj2(:,ColFlagIdx)=[];
    xDyT_BCAdj2(RowFlagIdx,:)=[];

    xDyT_BCAdj2=(xDyT_BCAdj2-min(xDyT_BCAdj2(:)))*1.0/(max(xDyT_BCAdj2(:))-min(xDyT_BCAdj2(:)));      
    xDyT_BCAdj2(isnan(xDyT_BCAdj2))=0;
    
    % imwrite(double(xDyT_BCAdj2),strcat(SavexDyTBCAdj0MeanFilePath,'xDyT_BCAdj2 -',num2str(k,'%04.0f'),'.tif'));

    xDyT_BCAdj0Mean=xDyT_BCAdj2-mean(xDyT_BCAdj2(:));
    
    HL_xDistyTime_BCAdj0Mean_Cell{k,1}=xDyT_BCAdj0Mean;

end

clearvars NumHistBins HistCtEdge HistCt_BinLoc HistCt BkgdHistBinIdx;
clearvars xDyT xDyT_Temp* xDyT_Col* xDyT_Row* xDyT_BCAdj*;
clearvars xDyT_Bkgd xDyT_Bkgd_FlagIdx;
clearvars ColFlagIdx RowFlagIdx ii;
clearvars *_Input;
 
 
%% Radon Transform: Slope Estimate
 
fprintf('Slope Estimate: Radon Transform...\n');
 
xDyT_Cell_Input=HL_xDistyTime_BCAdj0Mean_Cell;
 
Radon_Theta=-90:0.25:89.75;
 
Radon_Stdev_Mtx=zeros(numel(Radon_Theta),size(xDyT_Cell_Input,1));
 
HL_xDistyTime_RadonSlopeTheta=zeros(size(xDyT_Cell_Input,1),1);
 
for k=1:size(xDyT_Cell_Input,1) 
    
    [Radon_R,Radon_x]=radon(xDyT_Cell_Input{k,1},Radon_Theta);
    
    Radon_Stdev=std(Radon_R,0,1); 
       
    [~,Radon_MaxTheta_Idx]=max(Radon_Stdev(:));
    RadonSlopeTheta=Radon_Theta(Radon_MaxTheta_Idx)+90;
    
    Radon_Stdev_Mtx(:,k)=Radon_Stdev';
    HL_xDistyTime_RadonSlopeTheta(k,1)=RadonSlopeTheta;
    
end

clearvars HL_xDistyTime_BCAdj0Mean_Cell;
clearvars Radon_R Radon_x;
clearvars Radon_Stdev Radon_MaxTheta_Idx RadonSlopeTheta;
clearvars Radon_Theta;
clearvars *_Input;
 
 
%% Quality of Radon Slope Estimates
 
fprintf('Determine Quality of Radon Slope Estimate...\n');
 
Radon_Stdev_Mtx_Input=Radon_Stdev_Mtx;
 
[Radon_Stdev_Sort,Radon_Stdev_SortIdx]=sort(Radon_Stdev_Mtx_Input,1,'ascend');
Radon_Stdev_Bkgd=mean(Radon_Stdev_Sort(round(size(Radon_Stdev_Sort,1)*0.05):round(size(Radon_Stdev_Sort,1)*0.25),:),1);
 
Radon_Stdev2=Radon_Stdev_Mtx_Input./repmat(Radon_Stdev_Bkgd,size(Radon_Stdev_Mtx_Input,1),1);
 
Radon_Stdev2_Smooth=smoothdata(Radon_Stdev2,1,'gaussian',5);    % Smooth over 1 deg range
Radon_Stdev2_Max=max(Radon_Stdev2_Smooth,[],1);
 
HistCtEdge=1:0.25:ceil(max(Radon_Stdev2_Max(:)));
HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
[Radon_Stdev2_Max_HistCt,Radon_Stdev2_Max_HistCtEdge,Radon_Stdev2_Max_HistCtBin]=...
    histcounts(Radon_Stdev2_Max,HistCtEdge);
 
clearvars Radon_Stdev_Mtx;
clearvars Radon_Stdev_Sort Radon_Stdev_SortIdx;
clearvars Radon_Stdev_Bkgd;
clearvars Radon_Stdev2 Radon_Stdev2_Smooth;
clearvars HistCtEdge HistCt_BinLoc;
clearvars Radon_Stdev2_Max_HistCt Radon_Stdev2_Max_HistCtEdge;
clearvars *_Input;
 
 
%% Prepare Slope and Speed Data Matrix
 
fprintf('Prepare Data matrix...\n');
 
xDyT_Cell_Input=HL_xDistyTime_Cell;
xDyT_RadonSlopeTheta_Input=HL_xDistyTime_RadonSlopeTheta;
Radon_Stdev2_Max_Input=Radon_Stdev2_Max;
Radon_Stdev2_Max_HistCtBin_Input=Radon_Stdev2_Max_HistCtBin;

Flag_RadonSlopeTheta_0to90= find((0<=HL_xDistyTime_RadonSlopeTheta) & (HL_xDistyTime_RadonSlopeTheta<=90));
Flag_RadonSlopeTheta_90to180= find((90<HL_xDistyTime_RadonSlopeTheta) & (HL_xDistyTime_RadonSlopeTheta<=180));
 
HL_xDistyTime_RadonSlope1=atan(deg2rad(HL_xDistyTime_RadonSlopeTheta));
HL_xDistyTime_RadonSlope2=-atan(deg2rad(180-HL_xDistyTime_RadonSlopeTheta));
HL_xDistyTime_RadonSlope=zeros(size(HL_xDistyTime_RadonSlopeTheta,1),1);
HL_xDistyTime_RadonSlope(Flag_RadonSlopeTheta_0to90,1)=HL_xDistyTime_RadonSlope1(Flag_RadonSlopeTheta_0to90,1);
HL_xDistyTime_RadonSlope(Flag_RadonSlopeTheta_90to180,1)=HL_xDistyTime_RadonSlope2(Flag_RadonSlopeTheta_90to180,1);
clearvars HL_xDistyTime_RadonSlope1 HL_xDistyTime_RadonSlope2;
 
HL_xDistyTime_RadonSpeed=1./(HL_xDistyTime_RadonSlope.*Z_PxLength./XY_PxLength); 
 
% Build Data Matrix
% Each Row corresponds to 1 Hough Line
% Col1: Index of Hough Line (number assigned @ Hough Line generation)
% Col2: Radon Slope in Theta
% Col3: Radon Slope in y_px/x_px
% Col4: Radon Speed in um/s
% Col5: Radon Speed Est Quality, in Radon_Stdev2_Max
% Col6: Radon Speed Est Quality, in Binned and Ranked Radon_Stdev2_Max
xDyT_SlopeSpeed=horzcat((1:1:size(xDyT_Cell_Input,1))',...
    HL_xDistyTime_RadonSlopeTheta,...
    HL_xDistyTime_RadonSlope,...
    HL_xDistyTime_RadonSpeed,...
    Radon_Stdev2_Max',...
    Radon_Stdev2_Max_HistCtBin');

clearvars HL_xDistyTime_RadonSlopeTheta Radon_Stdev2_Max Radon_Stdev2_Max_HistCtBin;
clearvars Flag_RadonSlopeTheta_0to90 Flag_RadonSlopeTheta_90to180;
clearvars HL_xDistyTime_RadonSpeed;
clearvars HL_xDistyTime_RadonSlope HL_xDistyTime_RadonSlope1 HL_xDistyTime_RadonSlope2;
clearvars *_Input;
 
 
%% Overlay calculated slope on xDistyTime Images
 
fprintf('Overlap calculated xDist-yTimeDelay slope on xDistyTime images...\n');
 
xDyT_Cell_Input=HL_xDistyTime_Cell;
xDyT_SlopeSpeed_Input=xDyT_SlopeSpeed;
 
HL_xDistyTime_LinFitEst_Cell=cell(size(xDyT_Cell_Input,1),1);
 
for k=1:size(xDyT_Cell_Input,1)

    Img_xDyT=xDyT_Cell_Input{k,1}+abs(min(xDyT_Cell_Input{k,1}(:)));
    Img_xDyT=(Img_xDyT-min(Img_xDyT(:)))*0.75/(max(Img_xDyT(:))-min(Img_xDyT(:)));
 
    x=(1:1:size(xDyT_Cell_Input{k,1},2))';
    y=size(xDyT_Cell_Input{k,1},1)-round(xDyT_SlopeSpeed_Input(k,3).*x);
    y=y-min(y(:))+1;   
    x(y==0)=[];
    y(y==0)=[];
    x(y>size(Img_xDyT,1))=[];
    y(y>size(Img_xDyT,1))=[];
    xyLinIdx=sub2ind(size(Img_xDyT),y,x);
    Mask=zeros(size(Img_xDyT));
    Mask(xyLinIdx)=1;
    Img_xDyT_LinFitEst=max(Img_xDyT,Mask);
 
    HL_xDistyTime_LinFitEst_Cell{k,1}=Img_xDyT_LinFitEst;
    
end
 
clearvars HL_xDistyTime_Cell;
clearvars Img_xDyT x y xyLinIdx Mask Img_xDyT_LinFitEst;
clearvars *_Input;
 
 
%% Save HL_xDistyTime_LinFit images
 
fprintf('Saving xDistyTime_LinFit images...\n');
 
HL_xDistyTimeLinFitEst_Cell_Input=HL_xDistyTime_LinFitEst_Cell;
 
close all;
 
SavexDyTLFEFilePath=strcat(SaveFilePath,'xDistyTime_LinFitSlopeQual/');
mkdir(SavexDyTLFEFilePath);
 
for k=1:size(HL_xDistyTimeLinFitEst_Cell_Input,1)   
   imwrite(double(HL_xDistyTimeLinFitEst_Cell_Input{k,1}),strcat(SavexDyTLFEFilePath,'xDyTLinFit -',num2str(k,'%04.0f'),...
    ' S',num2str(xDyT_SlopeSpeed(k,4),'%0.2f'),...
    ' Q',num2str(xDyT_SlopeSpeed(k,6),'%0.2d'),...  
    '.tif'));
end

clearvars HL_xDistyTime_LinFitEst_Cell;
clearvars *_Input;
 
 
%% Map flow speed 
 
fprintf('Mapping flow speeds data...\n');
 
Img_PreSkel_Input=Img_PreSkel;  
 
HL_Mask_Stack_Input=HL_Mask_Stack;
xDyT_SlopeSpeed_Input=xDyT_SlopeSpeed;
 
VesselPx_LinIdx=find(Img_PreSkel_Input);
 
[~,xDyT_SlopeSpeed_Stdev2MaxSortIdx]=sortrows(xDyT_SlopeSpeed_Input,5,'descend');
HL_Rank=sortrows(horzcat(xDyT_SlopeSpeed_Stdev2MaxSortIdx,xDyT_SlopeSpeed_Input(:,1)),1);
 
ED_Stack=zeros(size(Img_PreSkel_Input,1),size(Img_PreSkel_Input,2),size(HL_Mask_Stack_Input,3));
 
for k=1:size(HL_Mask_Stack_Input,3)
        ED=bwdist(HL_Mask_Stack_Input(:,:,k));  
        ED_Stack(:,:,k)=ED;
end
 
EDDistLimit=zeros(size(HL_Mask_Stack_Input,3),1);
 
ED_Img_PreSkel_Neg=bwdist(~logical(Img_PreSkel_Input));
for k=1:size(HL_Mask_Stack_Input,3)
        EDDistLimit_Mtx=ED_Img_PreSkel_Neg.*logical(HL_Mask_Stack_Input(:,:,k));
        EDDistLimit_Mtx(EDDistLimit_Mtx==0)=NaN;
        EDDistLimit(k,1)=max(EDDistLimit_Mtx(:));    
end
 
EDDistLimit=repmat(min(EDDistLimit,[],1),size(EDDistLimit,1),1);

ED_Stack(~logical(Img_PreSkel_Input))=sqrt(size(Img_PreSkel_Input,1).^2+size(Img_PreSkel_Input,2).^2)+1;
ED_Stack(ED_Stack>repmat(permute(EDDistLimit,[2,3,1]),size(Img_PreSkel_Input,1),size(Img_PreSkel_Input,2),1))=...
    sqrt(size(Img_PreSkel_Input,1).^2+size(Img_PreSkel_Input,2).^2)+1;
 
Speed_HLIdx_Map=NaN(size(Img_PreSkel_Input));
Speed_Rank_Map=NaN(size(Img_PreSkel_Input));
Speed_Map=NaN(size(Img_PreSkel_Input));
 
Flag_HL_Signal1Rdm0=zeros(size(HL_Mask_Stack_Input,3),1);
Flag_HL_Signal1Rdm0(xDyT_SlopeSpeed_Input(:,6)>3)=1;
 
CellLength=5/XY_PxLength;
EDLVDistLimit=LargestVesselSize/XY_PxLength*0.5;
CellDistEdge=0:CellLength:EDLVDistLimit+CellLength;
CellDistEdge=horzcat(CellDistEdge,sqrt(size(Img_PreSkel_Input,1).^2+size(Img_PreSkel_Input,2).^2)+1);
 
for VesselPxCt=1:length(VesselPx_LinIdx)

    [RowIdx,ColIdx]=ind2sub([size(Img_PreSkel_Input,1),size(Img_PreSkel_Input,2)],VesselPx_LinIdx(VesselPxCt,1));
    
    Px_ED=permute(ED_Stack(RowIdx,ColIdx,:),[3,1,2]);  
   
    [Px_ED_HistCt,Px_ED_HistEdge,Px_ED_HistBin]=histcounts(Px_ED,CellDistEdge);
    
    Px_ED_HistBin(Px_ED_HistBin==length(CellDistEdge)-1)=9999;    
    Px_ED_EDDistLimit=ones(size(Px_ED_HistBin));   
    Px_ED_EDDistLimit(Px_ED_HistBin==9999)=0;
    
    [Sort,SortIdx]=sortrows(horzcat(Flag_HL_Signal1Rdm0,Px_ED_EDDistLimit,HL_Rank(:,2)),[1,2,3],{'descend','descend','ascend'});
 
     if Sort(1,1)==1 && Sort(1,2)>0 
        Speed_HLIdx_Map(VesselPx_LinIdx(VesselPxCt,1))=SortIdx(1,1);
        Speed_Rank_Map(VesselPx_LinIdx(VesselPxCt,1))=xDyT_SlopeSpeed_Input(SortIdx(1,1),6);      
        Speed_Map(VesselPx_LinIdx(VesselPxCt,1))=xDyT_SlopeSpeed_Input(SortIdx(1,1),4);           
     else       
        Speed_HLIdx_Map(VesselPx_LinIdx(VesselPxCt,1))=NaN;
        Speed_Rank_Map(VesselPx_LinIdx(VesselPxCt,1))=NaN;
        Speed_Map(VesselPx_LinIdx(VesselPxCt,1))=NaN;
     end
    
end
 
Speed_HLIdx_Map_flipud=flipud(Speed_HLIdx_Map);
Speed_Rank_Map_abs_flipud=abs(flipud(Speed_Rank_Map));
Speed_Map_abs_flipud=abs(flipud(Speed_Map));
 
% figHandle09=figure('Position',[10 150 1200 750]);
% hold on;
% axis on;
% axis equal;
% axis tight;
% plotHandle09=pcolor(Speed_HLIdx_Map_flipud);
% ylabel(strcat('Pixel position (y-axis)')); 
% xlabel('Pixel position (x-axis)');
% set(plotHandle09,'EdgeColor','none');
% colormap('lines');
% colorbar;
% ax=gca;
% ax.FontSize=18;
% title('Hough Line Index for Absolute Flow Speed Map','FontSize',12);
% hold off;
% print(figHandle09,'-dtiffn','-r0',strcat(SaveFilePath,'HoughLineldxMap r0.tif'));
 
figHandle11=figure('Position',[110 150 1200 750]);
hold on;
axis on;
axis equal;
axis tight;
plotHandle11=pcolor(Speed_Map_abs_flipud);
ylabel(strcat('Pixel position (y-axis)')); 
xlabel('Pixel position (x-axis)');
set(plotHandle11, 'EdgeColor' ,'none');
colormap('jet');
caxis([min(Speed_Map_abs_flipud(:)) 600]);
colorbar;
ax=gca;
ax.FontSize=18;
title('Absolute Flow Speed (um/s)','FontSize',12);
hold off;
print(figHandle11,'-dtiffn','-r300',strcat(SaveFilePath,'Abs FlowSpeed Map1 r300.tif'));
  
clearvars ImgStack;
clearvars HL_Mask_Stack;
clearvars *ED*;
clearvars xDyT_SlopeSpeed_Stdev2MaxSortIdx;
clearvars Flag_HL_Signal1Rdm0;
clearvars VesselPx_LinIdx;
clearvars HL_Rank;
clearvars Sort SortIdx;
clearvars CellLength CellDistEdge;
clearvars RowIdx ColIdx VesselPxCt k;
clearvars timestamp;
clearvars Save*Path -except SaveFilePath;
clearvars *_Input;
    
 
%% Save Flow Speed Data
 
fprintf('Saving xDyT_SlopeSpeed data...\n');
 
save(strcat(SaveFilePath,'/Speed_Map_abs_flipud.mat'),'Speed_Map_abs_flipud');
save(strcat(SaveFilePath,'/Speed_HLIdx_Map_flipud.mat'),'Speed_HLIdx_Map_flipud');
save(strcat(SaveFilePath,'/xDyT_SlopeSpeed.mat'),'xDyT_SlopeSpeed');
 
 
%% Save script in directory
 
% ScriptName=mfilename;
% PublishOptions=struct('format','html','showCode',true,'evalCode',false,'catchError',false,'figureSnapMethod','print','createThumbnail',false,'outputDir',SaveFilePath);
% publish(strcat(ScriptName,'.m'),PublishOptions);
 
 
%% 
% % % figure;
% % % for k=1:size(EDDistLimit,3)
% % %     imshow(double(EDDistLimit(:,:,k)));
% % %     % imshow(label2rgb(CellPos3D_LabelMtx(:,:,k)));
% % %     % imshow(Img_ZUnion);
% % %     % imwrite(double(ImgStack_Cav(:,:,k)),strcat(SaveCavFilePath,'Cav -',num2str(k,'%04.0f'),'.tif'));    
% % % end
