# This is plugin for ImageJ/Fiji. It is designed to assess the nuclear positivity of monitored antigens, both in H-DAB and fluorescence 
# stained samples. It is particularly useful for monitoring the translocation of nuclear receptors and other molecules from the cytoplasm 
# to the nucleus. The plugin evaluates only nuclear positivity and is not affected by cytoplasmic positivity.
# In addition to evaluation of nucleus-positive cells, the staining intensity in the nuclei can also be measured semiquantitatively 
# or quantitatively, depending on the staining method. 

# This script is for batch analysis enabling to analyze all images in folder. The results are written to .csv file and saved in the folder 
# with images (H-DAB) or in the folder of images with antibody signal (fluorescence / preprocessed images: H-DAB).
# Contrary to LoQANT for single analysis, only autothreshold for antibody signal is available.

# IMPORTANT NOTE: The correct view of .csv in Excell is dependent on separator used. This script use semicolon (;) as separator.
#				  Correct display of result can be set directly in Excell

# This program is distributed in the hope that it will be useful, but without any warranty.

# Author: Katerina Cizkova, July 2025


from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.gui import GenericDialog
from ij.plugin.frame import RoiManager 
from ij.measure import ResultsTable
from ij.io import DirectoryChooser
import os

# initial dialog for settings parameters:
def setting_parameters(): 
	gd = GenericDialog("Setting parameters for analysis")
	gd.addMessage("Choose relevant staining method: ")	
	staining = ["H-DAB", "fluorescence", "preprocessed images: H-DAB"]
	gd.addChoice("      Method: ", staining, staining[0])
	gd.addMessage("Image requierements:")
	gd.addMessage("* H-DAB: one RGB image")
	gd.addMessage("* Fluorescence: two RGB or grayscale images of separate channels")
	gd.addMessage("* Preprocessed images: H-DAB: two grayscale images of separate channels, white background")
	gd.addMessage("")
	gd.addMessage("Size of nuclei:")
	gd.addMessage("Set size of nuclei is highly recommended.")
	gd.addMessage("If there is no input, the range 0 - Infinity will be applied.")
	gd.addMessage("Set size of nuclei (px unit): ")
 	gd.addNumericField("         Minimum: ", 0, 0, 12, "")
 	gd.addNumericField("         Maximum: ", 1000000000, 0, 12, "")
	gd.addMessage("")
	intens = ["do not measure", "semiquantitative (rec for DAB)", "quantitative: positive nuclei (rec for fluorescence)","quantitative: all nuclei (rec for fluorescence)"]
	gd.addChoice("Nuclear intensity measurement:", intens, intens[0])
	gd.addMessage("")
	gd.addMessage("Additional settings:")
	segment_alg = ["Default", "Huang", "Huang2", "Intermodes", "IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"]
	gd.addChoice("Nuclear segmentation: ", segment_alg, segment_alg[0])
	gd.addChoice("Signal thresholding: ", segment_alg, segment_alg[0])
	perc = ["10", "20", "30", "40", "50", "60", "70", "80", "90", "100" ]
	gd.addChoice("                    Positive signal in nuclear area (%): ", perc, perc[7])
	gd.showDialog()
	method = gd.getNextChoice()
	minimum = gd.getNextNumber()
	maximum = gd.getNextNumber()
	if maximum == 1000000000:
		maximum = "Infinity"
	measure_intens = gd.getNextChoice()
	nucl_segm = gd.getNextChoice()
	signal_segm = gd.getNextChoice()
	overlap = float(gd.getNextChoice())/100
	if gd.wasCanceled():
		IJ.error("Analysis has been canceled by user.")
		next = False
	else:
		next = True
	return method, minimum, maximum, measure_intens, nucl_segm, signal_segm, overlap, next

def DAB_deconvolution(imp):
	IJ.run("Colour Deconvolution", "vectors=[H DAB] hide")
	imp_nuc = IJ.getImage(IJ.selectWindow(imp.title + "-(Colour_1)"))
	imp_sig = IJ.getImage(IJ.selectWindow(imp.title + "-(Colour_2)"))
	imp_rest = IJ.getImage(IJ.selectWindow(imp.title + "-(Colour_3)"))
	imp_rest.close()
	IJ.run(imp_nuc, "Set Scale...", "distance=1 known=0 unit=pixel")
	IJ.run(imp_sig, "Set Scale...", "distance=1 known=0 unit=pixel")
	return imp_nuc, imp_sig
	
def no_image():
	IJ.error("There is no image for analysis.")
	imp_nuc = None
	imp_sig = None
	next = False
	return imp_nuc, imp_sig, next
	
	
# detection of nuclei
def nuclei_detection(imp_nuclei, nucl_segment, minimum, maximum):
	IJ.run("Clear Results")
	IJ.run(imp_nuclei, "8-bit", "")
	if choosen_method == "H-DAB" or choosen_method == "preprocessed images: H-DAB":
		IJ.setAutoThreshold(imp_nuclei, "" + nucl_segment + "")
	else:
		IJ.setAutoThreshold(imp_nuclei, "" + nucl_segment + " dark")
	IJ.run(imp_nuclei, "Convert to Mask", "")
	IJ.run(imp_nuclei, "Watershed", "")
	IJ.run(imp_nuclei, "Analyze Particles...", "size=" + str(minimum) + "-" + str(maximum) + " exclude clear include add")
	imp_nuclei.show()
	IJ.run("Set Measurements...", "area mean limit display redirect=None decimal=3")
	global rm
	rm = RoiManager.getInstance()
	n = rm.getCount()
	for i in range(n):
		rm.getRoi(i)	
	rm.runCommand(imp_nuclei,"Measure") 
	table_nuclei = ResultsTable.getResultsTable()
	area_nuclei = []
	for i in range (n):
		area_nuclei += [table_nuclei.getValue("Area", i)]
	IJ.run("Close")
	return n, area_nuclei


# positive nuclear area
def measure_signal_area(imp_signal, n):
	IJ.run("Clear Results")
	IJ.run(imp_signal, "8-bit", "")
	if choosen_method == "H-DAB" or choosen_method == "preprocessed images: H-DAB":
		IJ.setAutoThreshold(imp_signal, "Default")
	else:
		IJ.setAutoThreshold(imp_signal, "Default dark")
	rm.runCommand(imp_signal,"Show All")
	rm.runCommand(imp_signal,"Measure")
	table_signal = ResultsTable.getResultsTable()
	area_signal=[]
	for i in range (n):
		area_signal += [table_signal.getValue("Area", i)]
	IJ.run("Close")                         
	return area_signal


def positive_nucl_area(area_signal, area_nuclei, n):
	positive_nuclear_area=[]
	for i in range (n):
		positive_nuclear_area += [area_signal[i]/area_nuclei[i]]
	return positive_nuclear_area


# percentage of positive cells
def positive_percentage(positive_nuclear_area, n, overlap):
	positive_nuclei_count = 0
	negative_nuclei_count = 0

	for i in range (n):
		if positive_nuclear_area[i] >= overlap:
			positive_nuclei_count += 1
		else:
			negative_nuclei_count += 1
	total = float(positive_nuclei_count + negative_nuclei_count)
	perc_positive_nuclei = (positive_nuclei_count*100)/(total)
	percentage_positive_nuclei = round((positive_nuclei_count*100)/(total), 3)
	return total, positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei


# quantitative intensity measurement
def quant_intens_positive(imp_signal, positive_nuclear_area, n, overlap):
	for i in range (n-1, -1, -1):
		if positive_nuclear_area[i] < overlap :
			rm.select(i)
			rm.runCommand(imp_signal, "Delete")
	reduced_n = rm.getCount()
	IJ.run("Clear Results")
	IJ.run("Set Measurements...", "area mean display redirect=None decimal=3")
	rm.runCommand(imp_signal,"Show All");
	rm.runCommand(imp_signal,"Measure");
	
	table_signal = ResultsTable.getResultsTable()
	intensity_signal=[]
	for i in range (reduced_n):
		mean_value = table_signal.getValue("Mean", i)
		intensity_signal += [round(mean_value, 3)]
	IJ.selectWindow("Results")
	IJ.run("Close")
	return intensity_signal

def quant_intens_all(imp_signal, n):
	IJ.run("Clear Results")
	IJ.run("Set Measurements...", "area mean display redirect=None decimal=3")
	rm.runCommand(imp_signal,"Show All")
	rm.runCommand(imp_signal,"Measure")
	table_signal = ResultsTable.getResultsTable()
	intensity_signal=[]
	for i in range (n):
		mean_value = table_signal.getValue("Mean", i)
		intensity_signal += [round(mean_value, 3)]
	IJ.run("Close")
	return intensity_signal


# DAB intensity measurement (semiquantitative)
def intensity_DAB(imp_signal, positive_nuclear_area, n, overlap):
	for i in range (n-1, -1, -1):
		if positive_nuclear_area[i] < overlap :
			rm.select(i)
			rm.runCommand(imp_signal, "Delete")
	reduced_n = rm.getCount()
	if reduced_n == 0:
		reduced_n = 0
		intensity_signal = []
	else:
		IJ.run("Clear Results")
		IJ.run("Set Measurements...", "area mean display redirect=None decimal=3")
		rm.runCommand(imp_signal,"Show All");
		rm.runCommand(imp_signal,"Measure");

		table_signal = ResultsTable.getResultsTable()
		intensity_signal=[]
		for i in range (reduced_n):
			intensity_signal += [table_signal.getValue("Mean", i)]
		IJ.selectWindow("Results")
		IJ.run("Close")
	return intensity_signal, reduced_n
	
def DAB_categories(positive_intensity_signal, n, reduced_n):	
	strong_positive_count = 0
	moderate_positive_count = 0
	weak_positive_count = 0
	total_count = float(n)
	total_DAB_intensity = ""
	
	if reduced_n == 0:
		perc_DAB_strong = 0
		perc_DAB_moderate = 0
		perc_DAB_weak = 0
		total_DAB_intensity = "Negative"
		hscore = 0
	else:	
		for i in range (reduced_n):
			if positive_intensity_signal[i] < 60:
				strong_positive_count += 1
			elif 60 <= positive_intensity_signal[i] < 120:
				moderate_positive_count += 1
			else: 
				weak_positive_count += 1
	
		total_count = float(n)
		perc_DAB_strong = round((strong_positive_count * 100)/total_count, 3)
		perc_DAB_moderate = round((moderate_positive_count * 100)/total_count, 3)
		perc_DAB_weak = round((weak_positive_count * 100)/total_count, 3)
		perc_DAB_neg = round(100-((reduced_n*100)/total_count), 3)
	
		total_score = perc_DAB_strong*4 + perc_DAB_moderate*3 + perc_DAB_weak*2 + perc_DAB_neg*1
	
		if perc_DAB_neg > 99:
			total_DAB_intensity = "Negative"
		elif total_score >= 301:
			total_DAB_intensity = "Strong positivity"
		elif (total_score >= 201) and (total_score < 301):
			total_DAB_intensity = "Moderate positivity"
		else:
			total_DAB_intensity = "Weak positivity"
		
		hscore = int(round(perc_DAB_strong*3 + perc_DAB_moderate*2 + perc_DAB_weak*1, 0))
		
	return strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, hscore 

# creating .csv files
def create_output_basic(folder_path, results):
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
	output_csv = os.path.join(folder_path, "LoQANT_basic_analysis_results.csv")
	with open(output_csv, "w") as f:
		for row in results:
			line = "{};{};{};{}; {}\n".format(row[0], row[1], row[2], row[3], row[4]) #space before {} is because of formating issue in Excell
			f.write(line)

				
def create_output_intensity(folder_path, measure_intensity, results):
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
	if measure_intensity == "semiquantitative (rec for DAB)":
		output_csv = os.path.join(folder_path, "LoQANT_semiquantitative_analysis_results.csv")
		with open(output_csv, "w") as f:
			for row in results:
				line = "{};{};{};{}; {};{};{};{}; {}; {}; {};{};{}\n".format(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12]) #spaces before {} are because of formating issue in Excell
				f.write(line)
				
	if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)":
		output_csv = os.path.join(folder_path, "LoQANT_quantitative_analysis_positive_results.csv")
		with open(output_csv, "w") as f:
			for row in results:
				line = "{}; {}\n".format(row[0], row[1]) #space before {} is because of formating issue in Excell
				f.write(line)

					
	if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
		output_csv = os.path.join(folder_path, "LoQANT_quantitative_analysis_all_results.csv")
		with open(output_csv, "w") as f:
			for row in results:
				line = "{}; {}\n".format(row[0], row[1]) #space before {} is because of formating issue in Excell
				f.write(line)

					

# close windows
def close_all(imp_nuclei, imp_signal):
	imp_nuclei.changes = False
	imp_nuclei.close()
	imp_signal.changes = False
	imp_signal.close()
	IJ.run("Clear Results")
	rm.close()


# basic analysis
choosen_method, min_size, max_size, measure_intensity, nuclei_segmentation, signal_segmentation, overlay, next_step = setting_parameters()
if next_step:
	if choosen_method == "H-DAB":
		dc = DirectoryChooser("Choose a folder with H-DAB images")
		folder_path = dc.getDirectory()

		if folder_path is None:
			IJ.error("The analysis has been canceled by user.")
			next_step = False
		
		if next_step:
			basic_results = [("image", "nuclei count", "positive", "negative", "% of positive")]
			semiquant_intensity_results = [("image", "nuclei count", "total positive", "negative", "% of positive", "strong positive", "moderate positive", "weak positive", "% strong positive", "% moderate positive", "% weak positive", "total DAB positivity", "histoscore")]
			quant_intensity_results = [("image", "mean gray value")]
				

			image_files = [f for f in os.listdir(folder_path) if f.lower().endswith((".tif", ".tiff", ".jpg", ".png"))]
			if image_files:
				for filename in image_files:
					full_path = os.path.join(folder_path, filename)
					imp = IJ.openImage(full_path)
					
					image_type = ImagePlus.getType(imp)
					if image_type == 4:	
						next_step = True
					else:
						imp_nuclei, imp_signal, next_step = no_image()
						
					if next_step:
						imp.show()
						imp_nuclei, imp_signal = DAB_deconvolution(imp)

						count_nuclei, area_nuclei = nuclei_detection(imp_nuclei, nuclei_segmentation, min_size, max_size)

						area_signal = measure_signal_area(imp_signal, count_nuclei)

						positive_nuclear_area = positive_nucl_area(area_signal, area_nuclei, count_nuclei)	

						total_count, positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei = positive_percentage(positive_nuclear_area, count_nuclei, overlay)
						
						basic_results.append((filename, int(total_count), positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei))
						
						if measure_intensity == "semiquantitative (rec for DAB)":
							positive_intensity_signal, reduced_n = intensity_DAB(imp_signal, positive_nuclear_area, count_nuclei, overlay)
							strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore = DAB_categories(positive_intensity_signal, count_nuclei, reduced_n)
							semiquant_intensity_results.append((filename, int(total_count), positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei, strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore))
						
						if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)":
							intensity_signal = quant_intens_positive(imp_signal, positive_nuclear_area, positive_nuclei_count, overlay)	
							for i in range (positive_nuclei_count):
								quant_intensity_results.append((filename, intensity_signal[i]))
						
						
						if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
							intensity_signal = quant_intens_all(imp_signal, count_nuclei)
							for i in range (count_nuclei):
								quant_intensity_results.append((filename, intensity_signal[i]))

						imp.close()
						close_all(imp_nuclei, imp_signal)						
						

			else:
				IJ.error("There are no images for analysis.")

			create_output_basic(folder_path, basic_results)
			if measure_intensity == "semiquantitative (rec for DAB)":
				create_output_intensity(folder_path, measure_intensity, semiquant_intensity_results)
			if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)": 
				create_output_intensity(folder_path, measure_intensity, quant_intensity_results)
			if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
				create_output_intensity(folder_path, measure_intensity, quant_intensity_results)


	else:
		dc_nuclei = DirectoryChooser("Choose a folder with images of nuclei")
		folder_path_nuclei = dc_nuclei.getDirectory()

		if folder_path_nuclei is None:
			IJ.error("The analysis has been canceled by user.")
			next_step = False
		
		dc_signal = DirectoryChooser("Choose a folder with images of signal")
		folder_path_signal = dc_signal.getDirectory()
		
		if folder_path_signal is None:
			IJ.error("The analysis has been canceled by user.")
			next_step = False
		
		if next_step:
			basic_results = [("image", "nuclei count", "positive", "negative", "% of positive")]
			semiquant_intensity_results = [("image", "nuclei count", "total positive", "negative", "% of positive", "strong positive", "moderate positive", "weak positive", "% strong positive", "% moderate positive", "% weak positive", "total DAB positivity", "histoscore")]
			quant_intensity_results = [("image", "mean gray value")]
							
			image_files_nuclei = [f for f in os.listdir(folder_path_nuclei) if f.lower().endswith((".tif", ".tiff", ".jpg", ".png"))]
			image_files_signal = [f for f in os.listdir(folder_path_signal) if f.lower().endswith((".tif", ".tiff", ".jpg", ".png"))]
			
			if image_files_nuclei and image_files_signal:

			# Loop through corresponding images from both folders
				for f1, f2 in zip(image_files_nuclei, image_files_signal):
					full_path_nuclei = os.path.join(folder_path_nuclei, f1)
					full_path_signal = os.path.join(folder_path_signal, f2)
					imp_nuclei = IJ.openImage(full_path_nuclei)
					imp_signal = IJ.openImage(full_path_signal)
					filename = imp_signal.getTitle()
					imp_nuclei.show()
					imp_signal.show()
					IJ.run(imp_nuclei, "Set Scale...", "distance=1 known=0 unit=pixel")
					IJ.run(imp_signal, "Set Scale...", "distance=1 known=0 unit=pixel")
					count_nuclei, area_nuclei = nuclei_detection(imp_nuclei, nuclei_segmentation, min_size, max_size)
					area_signal = measure_signal_area(imp_signal, count_nuclei)
					positive_nuclear_area = positive_nucl_area(area_signal, area_nuclei, count_nuclei)
					total_count, positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei = positive_percentage(positive_nuclear_area, count_nuclei, overlay)
					basic_results.append((filename, int(total_count), positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei))
					
					
					if measure_intensity == "semiquantitative (rec for DAB)" and choosen_method == "preprocessed images: H-DAB":
						positive_intensity_signal, reduced_n = intensity_DAB(imp_signal, positive_nuclear_area, count_nuclei, overlay)
						strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore = DAB_categories(positive_intensity_signal, count_nuclei, reduced_n)
						semiquant_intensity_results.append((filename, int(total_count), positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei, strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore))
					
					if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)":
						intensity_signal = quant_intens_positive(imp_signal, positive_nuclear_area, positive_nuclei_count, overlay)	
						for i in range (positive_nuclei_count):
							value = intensity_signal[i]
							rounded_value = round(float(value), 3)
							quant_intensity_results.append((filename, rounded_value))
					
					if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
						intensity_signal = quant_intens_all(imp_signal, count_nuclei)
						for i in range (count_nuclei):
							rounded_value = round(intensity_signal[i], 3)
							quant_intensity_results.append((filename, rounded_value))
					
					
					close_all(imp_nuclei, imp_signal)
					
				if measure_intensity == "semiquantitative (rec for DAB)" and choosen_method == "fluorescence":
					IJ.error("Semiquantitative analysis is not suitable. \n Only basic analysis will be performed.")
					next_step = False
					
		else:
			IJ.error("There are no images for analysis.")

		create_output_basic(folder_path_signal, basic_results)
		if next_step:
			if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)":
				create_output_intensity(folder_path_signal, measure_intensity, quant_intensity_results)
			if measure_intensity == "semiquantitative (rec for DAB)":
				create_output_intensity(folder_path_signal, measure_intensity, semiquant_intensity_results)
			if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
				create_output_intensity(folder_path_signal, measure_intensity, quant_intensity_results)
		

	

	