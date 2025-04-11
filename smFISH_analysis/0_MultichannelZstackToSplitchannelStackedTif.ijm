dir1 = getDirectory("Choose Source Directory"); 
list = getFileList(dir1);

for (i=0; i<list.length; i++ ) {
	thefile = dir1 + list[i];
	run("Bio-Formats Importer", "open=[thefile] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    filename = getTitle();
    filename_wo_ext = replace(filename, ".czi", ""); 
    foldername = filename_wo_ext + " TIFFs";
    dir2 = dir1 + File.separator + foldername;
    File.makeDirectory(dir2);
    
    run("Split Channels");

	stack1 = "C1-" + filename;
	selectImage(stack1);
	run("RGB Color");
	run("8-bit Color", "number=256");
	run("8-bit");
    saveAs("Tiff", dir1 + "/" + foldername + "/" + "C1-" + filename_wo_ext);
    
    stack2 = "C2-" + filename;
	selectImage(stack2);
    run("RGB Color");
    run("8-bit Color", "number=256");
	run("8-bit");
    saveAs("Tiff", dir1 + "/" + foldername + "/" + "C2-" + filename_wo_ext);
    
    stack3 = "C3-" + filename;
	selectImage(stack3);
    run("RGB Color");
    run("8-bit Color", "number=256");
	run("8-bit");
    saveAs("Tiff", dir1 + "/" + foldername + "/" + "C3-" + filename_wo_ext);
   
   run("Close All");
}