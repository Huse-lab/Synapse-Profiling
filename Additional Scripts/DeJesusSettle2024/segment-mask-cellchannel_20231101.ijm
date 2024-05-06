					   requires("1.33s"); 
					   dir = getDirectory("Choose a Directory ");
					   setBatchMode(true);
					   count = 0;
					   countFiles(dir);
					   n = 0;
					   processFiles(dir);
					   //print(count+" files processed");
					   
					   function countFiles(dir) {
					      list = getFileList(dir);
					      for (i=0; i<list.length; i++) {
					          if (endsWith(list[i], "/"))
					              countFiles(""+dir+list[i]);
					          else
					              count++;
					      }
					  }
					
					   function processFiles(dir) {
					      list = getFileList(dir);
					      for (i=0; i<list.length; i++) {
					          if (endsWith(list[i], "/"))
					              processFiles(""+dir+list[i]);
					          else {
					             showProgress(n++, count);
					             path = dir+list[i];
					             processFile(path);
					          }
					      }
					  }

// This is where you input the process to loop through all files.
function processFile(path) {
   if (matches(path, ".*.tif")) 
  	{
  		open(path);
		rename("orig");
		run("Split Channels");
		selectWindow("C1-orig");
		close();
		selectWindow("C2-orig");
		Stack.getDimensions(width, height, channels, slices, frames);
		// Perform by time point
		for (timepoint=1; timepoint<=frames; timepoint++){
			Stack.setFrame(timepoint);
			run("Duplicate...", "duplicate frames=" + timepoint);
			// Run operations
			run("Gaussian Blur 3D...", "x=2 y=2 z=2");
			setAutoThreshold("Default dark");
			//run("Threshold...");
			setThreshold(1, 65535);
			setOption("BlackBackground", true);
			run("Convert to Mask", "method=Default background=Dark black");
			run("Fill Holes", "stack");
			run("Gaussian Blur 3D...", "x=2 y=2 z=2");
			run("3D Objects Counter", "threshold=1 slice=1 min.=100000 max.=9999999999999 objects");
			filename_segmentmask = replace(path,".tif","_SegmentMask_");
			savename = filename_segmentmask + "time" + timepoint + ".tif";
			saveAs("Tiff",savename);
			run("Close");
			run("Close");
			}
		run("Close All");	
  	}
}
  
print("Done creating segment masks.");