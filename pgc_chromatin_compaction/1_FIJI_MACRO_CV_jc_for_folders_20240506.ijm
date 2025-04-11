dir0 = getDirectory("Choose main directory with subfolders whose names have Copy in the end"); 
regexstring = ".*" + "- Copy" + ".*"
list = getFileList(dir0);

for (m=0; m<list.length; m++) {
	if (matches(list[m], regexstring)) {
		// dir1 = dir0 + "/" + list[m] + "/";
		dir1 = dir0 + File.separator + list[m];
		imagelist = getFileList(dir1);
		
		for (i=0; i<imagelist.length; i++) {
			open(dir1 + imagelist[i]);
			imageName = getTitle();
			run("Enhance Contrast", "saturated=0.35");
			run("Duplicate...", " ");
			run("Gaussian Blur...", "sigma=2");
			setAutoThreshold("Default dark");
			setThreshold(4, 255); // JC changed lower limit. Always set upper limit to 255
			setOption("BlackBackground", false);
			run("Convert to Mask");
			run("Fill Holes");
			run("Dilate");
			run("Watershed");
			run("Analyze Particles...", "size=2000-Infinity pixel exclude add"); //JC changed size= and added pixel
			selectWindow(imageName);
			n = roiManager("count");
			
			if (n > 0) {
				for (j=0; j < n; j++) {
				roiManager("select", j);
				// title = getTitle();
				run("Duplicate...", " ");
				setBackgroundColor(0, 0, 0);
				run("Clear Outside");
				saveAs(dir1+imageName+"_"+(j+1)+"_ori.tif");
				
				run("Create Mask");
				saveAs(dir1+imageName+"_"+(j+1)+"_fill.tif"); // gray scale single-channel so black is 255 and white is 0
				// have to close twice because now we have an ROI window and a mask window
				close();				
				close();
				}
				
				selectWindow(imageName);
				run("Enhance Contrast", "saturated=0.35");
				roiManager("Show All with labels");
				// roiManager("Save", dir1+imageName+"_"+"ROI"+".zip");
				wait(1000);
				run("Capture Screen");
				saveAs(dir1+imageName+"_"+"Screenshot"+".jpg");
				close();
				run("Close All");
				
				selectWindow("ROI Manager");
				run("Close");
			} else {
				selectWindow("ROI Manager");
				run("Close");
				run("Close All");
			}			
		}
	}
}
	
