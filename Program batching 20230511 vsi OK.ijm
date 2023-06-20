ImgFolder=getDirectory("Import an Img folder");
Imglist=getFileList(ImgFolder);
SavingFolder=getDirectory("Saving folder");
CellSize=0;
//FitExt=newArray("tif","tiff");
FitExt=newArray("vsi");
for (i = 0; i < lengthOf(Imglist); i++){
	print(Imglist[i]);
	close("*");
	if(FitExtension(Imglist[i],FitExt)){
		ImgSavingFolder=SavingFolder+File.getNameWithoutExtension(Imglist[i])+"\\";
		File.makeDirectory(ImgSavingFolder);
		run("Bio-Formats", "open=["+ImgFolder+Imglist[i]+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		Img=getTitle();
		if(CellSize==0){
			arr=nChannelsAsArray(Img);
			Dialog.createNonBlocking("Channel selection");
			Dialog.addChoice("Microtubule", arr);
			Dialog.addChoice("Nuclear", arr);
			Dialog.addChoice("Mitochondria", arr);
			Dialog.addNumber("Cell size", 800);
			Dialog.show();
			Tub =Dialog.getChoice();
			Nu  =Dialog.getChoice();
			Mito=Dialog.getChoice();
			CellSize=Dialog.getNumber();		
		}
		selectWindow(Masking(Img,Nu,Mito,Tub));
		saveAs("Tiff", ImgSavingFolder+"Mask.tiff");
		MaskImg=getTitle();
		setSlice(1);
		run("Analyze Particles...", "size=[&CellSize-Infinity] pixel exclude clear include add slice");
		roiManager("save", ImgSavingFolder+"RoiSet.zip");
		File.makeDirectory(ImgSavingFolder+"ROIs\\");
		//Open nuclear image
		close("*");
		run("Bio-Formats", "open=["+ImgFolder+Imglist[i]+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		Channels=nChannelsAsArray(getTitle());
		Img=getTitle();
		run("Split Channels");
		for (ch = 0; ch <lengthOf(Channels); ch++) {
			if(Channels[ch]==Nu)NucImg="C"+Channels[ch]+"-"+Img;
			else close("C"+Channels[ch]+"-"+Img);
		}

		CellTrack(ImgSavingFolder+"Mask.tiff", NucImg, ImgSavingFolder+"ROIs\\", ImgSavingFolder+"RoiSet.zip", 2,CellSize);
		close("*");
		run("Bio-Formats", "open=["+ImgFolder+Imglist[i]+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		Img=getTitle();
		File.makeDirectory(ImgSavingFolder+"CSVs\\");
		getCSV(Img,ImgSavingFolder+"ROIs\\",ImgSavingFolder+"CSVs\\");
		//saveAs("tiff", ImgSavingFolder+Imglist[i]);
		//File.rename(ImgFolder+Imglist[i],ImgSavingFolder+Imglist[i]);
	}
}
print("\\Clear");
print("Nu,"+Nu);
print("Tub,"+Tub);
print("Mito,"+Mito);
print("Cellsize,"+CellSize);
saveAs("text",SavingFolder+"Log.csv");
close("Log");
close("*");

//--------------------------------------------
//Mask Function
function Masking(Img,Nu,Mito,Tub){
	selectWindow(Img);
	getDimensions(width, height, channels, slices, frames);
	run("Split Channels");
	for (i = 1; i <= channels; i++) {
		if(i!=Nu&&i!=Mito&&i!=Tub)close("C"+i+"-"+Img);
	}
	Tub ="C"+Tub+"-"+Img;
	Nu  ="C"+Nu+"-"+Img;
	Mito="C"+Mito+"-"+Img;
	thr1=BackgroundRemove(Tub,100);
	for(i=0;i<2;i++){
		run("Gaussian Blur...", "sigma=2 stack");
		run("Subtract Background...", "rolling=50 sliding stack");
	}
	thr1=threshold(thr1,0,true);		//Microtubular mask
	selectWindow(Mito);
	run("Subtract Background...", "rolling=50 sliding stack");
	thr2=threshold(Mito, 0,true);		//Mitochondrial mask
	thr1=Comb(thr1,thr2,"OR");			//Combine microtubular and mitochondrial mask
	selectWindow(Nu);
	run("Subtract Background...", "rolling=50 sliding stack");
	run("Duplicate...", "duplicate");
	rename("Temp");						//Nuclear Mask
	thr2=threshold("Temp",0,true);
	thr1=Comb(thr1,thr2,"OR");
	run("Fill Holes", "stack");
	close("Temp");
	selectWindow(Nu);
	seg=segment(Nu, 10, 50);			//Nuclear Segmentation
	thr1=Comb(thr1,seg,"AND");			//Segment combine
	selectWindow(thr1);
	close("\\Others");
	rename("Mask");
	return "Mask";
}
function Comb(Img1,Img2,Mod){
	Mod=toUpperCase(Mod);
	if(Mod=="AND")imageCalculator("AND stack",Img1,Img2);
	else if(Mod=="OR")imageCalculator("OR stack",Img1,Img2);
	else if(Mod=="DIV32"){
		imageCalculator("Divide creat 32-bit stack", Img1, Img2);
		close(Img1);
		rename(Img1);
	}
	close(Img2);
	return Img1;
}
function threshold(Window, thr, yAutoContrast){
	selectWindow(Window);
	if(yAutoContrast){
		selectWindow(Window);
		run("Enhance Contrast", "saturated=0.35");
		run("Apply LUT", "stack");
	}
	if(thr>0)setThreshold(thr,pow(2,bitDepth())-1);
	else setAutoThreshold("Default dark");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Dark");
	run("Fill Holes", "stack");
	run("Erode", "stack");
	for(i=0;i<3;i++)run("Dilate", "stack");
	run("Fill Holes", "stack");
	rename(Window+"_thr");
	return Window+"_thr";
}
function segment(Window,sigma,prominence){
	selectWindow(Window);
	n=nSlices;
	run("Make Montage...", "columns=1 rows=&n scale=1");
	selectWindow("Montage");
	run("Gaussian Blur...", "sigma=&sigma");
	run("Find Maxima...", "prominence=&prominence exclude output=[Segmented Particles]");
	close("Montage");
	selectWindow("Montage Segmented");
	run("Montage to Stack...", "columns=1 rows=&n border=0");
	close("Montage Segmented");
	selectWindow("Stack");
	rename("seg");
	return "seg";
}
function BackgroundRemove(Img,sigma){
	selectWindow(Img);
	run("Duplicate...", "title=[Back] duplicate");
	run("Gaussian Blur...","sigma=[&sigma] stack");
	Img=Comb(Img,"Back","DIV32");
	run("16-bit");
	rename(Img+"_BKRemove");
	return Img+"_BKRemove";
}

//Function Tracking
function CellTrack(MaskImg,NucImg,SVpath,RefROI,Radius,CellSize){//MaskImg: filepath, NucImg: Opening window name
	roiManager("reset");
	roiManager("Open", RefROI);
	RefCount=roiManager("count");
	for(i=0;i<RefCount;i++){
		roiManager("reset");
		roiManager("Open", RefROI);
		UnprocessingROIs= getDeletingROI(i);
		roiManager("select",UnprocessingROIs);
		roiManager("Delete");
		open(MaskImg);
		MaskWindow=getTitle();
		run("Properties...", "channels=1 slices=1 frames="+nSlices+" pixel_width=1 pixel_height=1 voxel_depth=1");
		TrackingCorr=Tracker(MaskWindow,NucImg,CellSize,Radius);
		RepeatROIs=getRepeatROIs(true);
		if(lengthOf(RepeatROIs)>0){
			roiManager("select", RepeatROIs);
			roiManager("delete");
		}
		close(MaskWindow);
		roiManager("save", SVpath+i+".zip");
		roiManager("reset");
	}
	close("Results");
	close("*");
	return MaskImg;
}
function Tracker(Mask,Raw,Size,Radius){
	selectWindow(Mask);
	Target=1;
	for(i=1;i<nSlices;i++){
		if(Target>0){
			roiManager("Select",i-1);
			r=Radius*2*(getValue("selection.size")/3.14159)^0.5;	//Analyze particles範圍
			x=getValue("X");
			y=getValue("Y");
			Target=ROITracker(Mask,i-1,i,Size,x,y,r,"size=[&Size-Inf] pixel display add slice");
			run("Select None");
		}
	}
	if(Target>0){
		MitoticFrame=getMitosis(Raw);
		if(MitoticFrame>0){
			setForegroundColor(255, 255, 255);
			selectWindow(Mask);
			roiManager("select", MitoticFrame);
			X_DaughterA=getValue("X");
			Y_DaughterA=getValue("Y");
			roiManager("Fill");
			roiManager("deselect");
			roiManager("select", MitoticFrame-1);
			X_Mother=getValue("X");
			Y_Mother=getValue("Y");
			x=X_Mother*2-X_DaughterA;
			y=Y_Mother*2-Y_DaughterA;
			r=Radius*2*(getValue("selection.size")/3.14159)^0.5;
			ROITracker(Mask,MitoticFrame-1,MitoticFrame,200,x,y,r,"size=[&Size-Inf] pixel display add slice exclude");
			for (i = MitoticFrame+1; i < nSlices; i++){
				roiManager("Select",roiManager("count")-1);
				r=Radius*2*(getValue("selection.size")/3.14159)^0.5;	//Analyze particles範圍
				x=getValue("X");
				y=getValue("Y");
				ROITracker(Mask,roiManager("count")-1,i,200,x,y,r,"size=[&Size-Inf] pixel display add slice");
				run("Select None");
			}
		}
	}
	return Target;
}
function ROITracker(Img,RefROI,Frame,Size,x,y,r,Command_str){
	selectWindow(Img);
	run("Select None");
	setSlice(Frame+1);
	makeOval(x-r/2,y-r/2, r, r);
	roiMan_n = roiManager("count");
	run("Analyze Particles...", Command_str);
	roiMan_n_new = roiManager("count");
	//計算下一張Slice增加了多少個ROI
	roi_add = roiMan_n_new - roiMan_n ;

	roiManager("Show None");
	AND_Big_area= 0; 	//預設取交集後的面積為零
	AND_Big_area_n = 0 ;//預設取交集後的面積，最大的是第零個
	if(roi_add>1){
		for (j=1; j<=roi_add; j++) {
			roiManager("Select", newArray(RefROI,roiMan_n-1+j));
			roiManager("AND");
			area=getValue("selection.size");
			//沒交集的面積會回報成整張圖片的面積，所以將回報面積設為0
			if(area == getWidth() * getHeight())area=0;
			else if(area>AND_Big_area){
				AND_Big_area = area;
				AND_Big_area_n = j;
			}
		}
		//select biggest area and add the new ROI in the next slice
		roiManager("Select", roiMan_n -1 + AND_Big_area_n);
		roiManager("Add");
	
		//delete all new ROI in the next slice
		for (j=1; j<=roi_add; j++) {
			roiManager("deselect");
			roiManager("Select", roiMan_n);
			roiManager("Delete");
		}
	}
	run("Clear Results");
	return roi_add;
}
function getMitosis(Img){
	run("Clear Results");
	selectWindow(Img);
	roiManager("deselect");
	roiManager("Measure");
	Arr=Table.getColumn("IntDen");
	run("Clear Results");
	MaxDeltaFrame=0;
	MaxDelta=0;
	for (i = 0; i < lengthOf(Arr)-1; i++) {
		Delta=Arr[i+1]-Arr[i];
		if(Delta<MaxDelta){
			//add Average Circ. 3 frames before mitosis > 0.75
			MaxDeltaFrame=i;
			MaxDelta=Delta;
		}
	}
	if(MaxDeltaFrame>0)return MaxDeltaFrame+1;
	else return 0;
}

function getDeletingROI(RemainROI){
	return Array.deleteIndex(Array.getSequence(roiManager("count")), RemainROI);
}
function getRepeatROIs(isStack){
	RepeatArr=newArray(0);
	run("Select None");
	for (i = 0; i < roiManager("count")-1; i++) {
		roiManager("select", i);
		Roi.getPosition(channel, slice, frame);
		slice_i=slice;
		for (j = i+1; j < roiManager("count"); j++) {
			roiManager("select", newArray(j,i));
			Roi.getPosition(channel, slice, frame);
			slice_j=slice;
			roiManager("XOR");
			getStatistics(area, mean, min, max, std, histogram);
			if(area==getWidth()*getHeight()){
				if(slice_j==slice_i||!isStack)RepeatArr=Array.concat(RepeatArr,j);
			}
		}
	}
	return RepeatArr;
}

//ROI to CSV Function
function getCSV(Img,ROIdir,SavePath){
	selectWindow(Img);
	getDimensions(width, height, channels, slices, frames);
	for (i = 1; i <= channels; i++)getCSVCh(Img,ROIdir,SavePath,i);
}
function getCSVCh(Img,ROIdir,SavePath,Ch){
	selectWindow(Img);
	Stack.setChannel(Ch);
	list=getFileList(ROIdir);
	File.makeDirectory(SavePath+"Ch"+Ch+"\\");
	SVPath=SavePath+"Ch"+Ch+"\\";
	for (i=0;i<list.length;i++) {
		roiManager("open",ROIdir+list[i]);
		roiManager("Measure");
		saveAs("Results", SVPath+"Result "+SerialNum(getID(list[i]), digSize(list.length))+".csv");
		close("Results");
		roiManager("reset");
	}
}
function getID(Name){
	Name=split(Name,".");
	return Name[0];
}
function SerialNum(i, dig){
	Temp="";
	Temp=Temp+i;
	if(lengthOf(Temp)<dig){
		while(lengthOf(Temp)<dig)Temp="0"+Temp;
	}
	return Temp;
}
function digSize(i){
	Temp="";
	Temp=Temp+i;
	return lengthOf(Temp);
}

//Function bank
function nChannelsAsArray(Img){
	selectWindow(Img);
	getDimensions(width, height, channels, slices, frames);
	arr="";
	for (i = 1; i <=channels; i++)arr=arr+" "+i;
	return split(arr, " ");
}
function FitExtension(Filename,fitExt){
	FileExt=getExtension(Filename);
	for (i = 0; i < lengthOf(fitExt); i++) {
		if(matches(FileExt, fitExt[i]))return true;
	}
	return false;
}
function getExtension(Path){
	Arr= split(Path,".");
	return toLowerCase(Arr[lengthOf(Arr)-1]);
}
