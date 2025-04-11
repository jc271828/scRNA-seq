dir1 = getDirectory("Choose Source Directory"); 
list = getFileList(dir1);
for (i=0; i<list.length; i++ ) {
	thefile = dir1 + list[i];
	run("Bio-Formats Importer", "open=[thefile] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    //run("Bio-Formats Importer", "open=" + dir1 + what + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    //open(dir1 + list[i]);
    filename = getTitle();
    filename_wo_ext = replace(filename, ".czi", ""); 
    foldername = filename_wo_ext + " TIFFs";
    dir2 = dir1 + File.separator + foldername;
    File.makeDirectory(dir2);
        
    run("Hyperstack to Stack");
    run("RGB Color");
    run("8-bit Color", "number=256");
    run("Stack to Images");
    
    ids = newArray(nImages);
    for (j=0; j<nImages; j++) { 
        selectImage(j + 1); 
        title = getTitle();
        //saveAs("Tiff", dir2 + title); somehow directly using dir2 doesn't work
        saveAs("Tiff", dir1 + "/" + foldername + "/" + title);
        }
    
    run("Close All");
}