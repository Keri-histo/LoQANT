# This is plugin for ImageJ/Fiji. It is designed to assess the nuclear positivity of monitored antigens, both in H-DAB and fluorescence 
# stained samples. It is particularly useful for monitoring the translocation of nuclear receptors and other molecules from the cytoplasm 
# to the nucleus. The plugin evaluates only nuclear positivity and is not affected by cytoplasmic positivity.
# In addition to evaluation of nucleus-positive cells, the staining intensity in the nuclei can also be measured semiquantitatively 
# or quantitatively, depending on the staining method. 

# This program is distributed in the hope that it will be useful, but without any warranty

# Author: Katerina Cizkova, January 2025

from ij import IJ, ImagePlus
from ij import WindowManager as WM
from ij.gui import GenericDialog
from ij.plugin.frame import RoiManager 
from ij.measure import ResultsTable
from java.awt import event, Font
from java.awt.event import AdjustmentListener 
from javax.swing import JFrame, JTextArea
import os

# staining method
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
	
# images for analysis
def select_images_for_analysis(method):
	if method == "H-DAB":
		if WM.getImageCount() == 0:
			imp_nuc, imp_sig, next = no_image()
		else:
			imp = IJ.getImage()
			image_type = ImagePlus.getType(imp)
			if image_type == 4:
				imp_nuc, imp_sig = DAB_deconvolution(imp)
				next = True
			else:
				imp_nuc, imp_sig, next = no_image()
	else:
		if WM.getImageCount() < 2:
			imp_nuc, imp_sig, next = no_image()
		else:
			imp = IJ.getImage()
			image_titles_unicode = WM.getImageTitles()
			image_titles = [title.encode('utf-8') for title in image_titles_unicode]
			image_options = process_images_titles(image_titles)
			imp_nuc, imp_sig, next = select_channels(image_options)
			IJ.run(imp_nuc, "Set Scale...", "distance=1 known=0 unit=pixel")
			IJ.run(imp_sig, "Set Scale...", "distance=1 known=0 unit=pixel")
	return imp_nuc, imp_sig, next

# fluorescent or preprocessed images channel selection
def process_images_titles(titles):
	image_options = []
	for title in titles:
		image_options.append(title)
	return image_options

def select_channels(image_options): 
	gd = GenericDialog("Choose channels")
	gd.addMessage("Images of separate channels are needed.")
	gd.addMessage("Choose the channels properly:")
	gd.addChoice("nuclei: ", image_options, image_options[0])
	gd.addChoice("signal: ", image_options, image_options[0])
	gd.showDialog()
	nuclei = gd.getNextChoice()
	signal = gd.getNextChoice()
	imp_nucl = IJ.getImage(IJ.selectWindow(nuclei))
	IJ.run(imp_nucl, "8-bit", "")
	imp_sig = IJ.getImage(IJ.selectWindow(signal))
	IJ.run(imp_sig, "8-bit", "")
	if gd.wasCanceled():
		IJ.error("Analysis has been canceled by user.")
		next = False
	else:
		next = True
	return imp_nucl, imp_sig, next


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


# antibody signal threshold
def threshold_image_UI(auto_value, method):
	gd = GenericDialog("Threshold for antibody signal")
	gd.addMessage("The threshold has been set automatically to: " + str(auto_value))
	gd.addMessage("Use of this settings is highly recommended.")
	gd.addMessage("")
	gd.addMessage("If it is necessary, the threshold can be set manually:")
	gd.addSlider("   ", 0, 255, auto_value)
	slider = gd.getSliders()[0]
	if method == "H-DAB" or method == "preprocessed images: H-DAB":
		slider.addAdjustmentListener(lambda event: threshold_changed_DAB(slider))
	else:
		slider.addAdjustmentListener(lambda event: threshold_changed_fluorescence(slider))	
	gd.addMessage("")
	gd.showDialog()
	if gd.wasCanceled():
		IJ.error("Analysis has been canceled by user.")
		next = False
		rm.close()
	else:
		next = True
	return next, slider

def threshold_changed_DAB(slider):
    threshold_value = slider.getValue() 
    IJ.setThreshold(imp_signal, 0, threshold_value) 
    imp_signal.updateAndDraw()	
    return threshold_value

def threshold_changed_fluorescence(slider):
    threshold_value = slider.getValue()
    IJ.setThreshold(imp_signal, threshold_value, 255)
    imp_signal.updateAndDraw()
    return threshold_value

def signal_threshold(imp_signal):
	IJ.run(imp_signal, "8-bit", "")
	if choosen_method == "H-DAB" or choosen_method == "preprocessed images: H-DAB":
		IJ.setAutoThreshold(imp_signal, "Default")
	else:
		IJ.setAutoThreshold(imp_signal, "Default dark")
	min_threshold = imp_signal.getProcessor().getMinThreshold()
	max_threshold = imp_signal.getProcessor().getMaxThreshold()
	if choosen_method == "H-DAB" or choosen_method == "preprocessed images: H-DAB":
		next, slider = threshold_image_UI(max_threshold, choosen_method)
	else:
		next, slider = threshold_image_UI(min_threshold, choosen_method)
	return next


# positive nuclear area
def measure_signal_area(imp_signal, n):
	IJ.run("Clear Results")
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
	percentage_positive_nuclei = round((positive_nuclei_count*100)/(total), 2)
	return total, positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei


def show_results(total, negative, positive, percentage):
	frame = JFrame("Results of counting")
	text_area = JTextArea("Total cell count: " +str(int(total))+
	           "\nNegative nuclei: " +str(negative)+ 
	           "\nPositive nuclei: " +str(positive)+ 
	           "\n\nPercentage of positive nuclei: " +str(percentage)+ " %")
	text_area.font = Font("Calibri", Font.PLAIN, 15)
	frame.add(text_area)
	frame.setSize(300,200)
	frame.setLocation(250,300)
	frame.setVisible(True)


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
		intensity_signal += [table_signal.getValue("Mean", i)]
	IJ.selectWindow("Results")
	IJ.run("Close")
	return intensity_signal

def quant_intens_all(imp_signal, n):
	IJ.run("Clear Results")
	IJ.run("Set Measurements...", "area mean display redirect=None decimal=3")
	rm.runCommand(imp_signal,"Show All")
	rm.runCommand(imp_signal,"Measure")
	table_signal = ResultsTable.getResultsTable()
	intens_signal=[]
	for i in range (n):
		intens_signal += [table_signal.getValue("Mean", i)]
	IJ.run("Close")                         
	return intens_signal

def show_intensity_tab_all(signal_intensity):
	rt_int_all = ResultsTable()
	rt_int_all.setHeading(0, "No.")
	rt_int_all.setHeading(1, "mean gray value")
	label = 1
	for intensity in signal_intensity:
		rt_int_all.incrementCounter()
		rt_int_all.addValue("No.", label)
		rt_int_all.addValue("mean gray value", intensity)
		label += 1
	rt_int_all.show("Intensity of all nuclei")
	
def show_intensity_tab_positive(signal_intensity):
	rt_int_positive = ResultsTable()
	rt_int_positive.setHeading(0, "No.")
	rt_int_positive.setHeading(1, "mean gray value")
	label = 1
	for intensity in signal_intensity:
		rt_int_positive.incrementCounter()
		rt_int_positive.addValue("No.", label)
		rt_int_positive.addValue("mean gray value", intensity)
		label += 1
	rt_int_positive.show("Intensity of positive nuclei")

# DAB intensity measurement
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
		IJ.run("Set Measurements...", "area mean limit display redirect=None decimal=3")
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
		perc_DAB_strong = round((strong_positive_count * 100)/total_count, 1)
		perc_DAB_moderate = round((moderate_positive_count * 100)/total_count, 1)
		perc_DAB_weak = round((weak_positive_count * 100)/total_count, 1)
		perc_DAB_neg = round(100-((reduced_n*100)/total_count), 1)
	
		total_score = perc_DAB_strong*4 + perc_DAB_moderate*3 + perc_DAB_weak*2 + perc_DAB_neg*1
	
		if perc_DAB_neg > 99:
			total_DAB_intensity = "Negative"
		elif total_score >= 301:
			total_DAB_intensity = "Strong positivity"
		elif (total_score >= 201) and (total_score < 301):
			total_DAB_intensity = "Moderate positivity"
		else:
			total_DAB_intensity = "weak positivity"
		
		hscore = int(round(perc_DAB_strong*3 + perc_DAB_moderate*2 + perc_DAB_weak*1, 0))
		
	return strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, hscore 

def show_results_DAB_intensity(positive, DAB_strong, DAB_moderate, DAB_weak, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, hscore):
	frame = JFrame("Results of DAB intensity")
	text_area = JTextArea("Positive nuclei: " +str(positive)+ "\nStrong positivity: " +str(DAB_strong)+ "\nModerate positivity: " +str(DAB_moderate)+ "\nWeak positivity: " +str(DAB_weak)+ 
		"\n\nPercentage of strong positivity: " +str(perc_DAB_strong)+ 
		" %\nPercentage of moderate positivity: " +str(perc_DAB_moderate)+ " %\nPercentage of weak positivity: " +str(perc_DAB_weak)+ 
		" %\n\nTotal positivity: " +str(total_DAB_intensity)+ " \n\nHistoscore: " +str(hscore))
	text_area.font = Font("Calibri", Font.PLAIN, 15)
	frame.add(text_area)
	frame.setSize(300,300)
	frame.setLocation(600, 300)
	frame.setVisible(True)



# basic analysis
choosen_method, min_size, max_size, measure_intensity, nuclei_segmentation, signal_segmentation, overlay, next_step = setting_parameters()
if next_step:
	imp_nuclei, imp_signal, next_step = select_images_for_analysis(choosen_method)
	
	if next_step:
		count_nuclei, area_nuclei = nuclei_detection(imp_nuclei, nuclei_segmentation, min_size, max_size)
		next_step = signal_threshold(imp_signal)
		
		if next_step:
			area_signal = measure_signal_area(imp_signal, count_nuclei)

			positive_nuclear_area = positive_nucl_area(area_signal, area_nuclei, count_nuclei)	

			total_count, positive_nuclei_count, negative_nuclei_count, percentage_positive_nuclei = positive_percentage(positive_nuclear_area, count_nuclei, overlay)

			show_results(total_count, negative_nuclei_count, positive_nuclei_count, percentage_positive_nuclei)


# intensity measurement
if measure_intensity == "semiquantitative (rec for DAB)":
	positive_intensity_signal, reduced_n = intensity_DAB(imp_signal, positive_nuclear_area, count_nuclei, overlay)
	strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore = DAB_categories(positive_intensity_signal, count_nuclei, reduced_n)
	show_results_DAB_intensity(positive_nuclei_count, strong_positive_count, moderate_positive_count, weak_positive_count, perc_DAB_strong, perc_DAB_moderate, perc_DAB_weak, total_DAB_intensity, histoscore)

if measure_intensity == "quantitative: all nuclei (rec for fluorescence)":
	intensity_signal = quant_intens_all(imp_signal, count_nuclei)
	show_intensity_tab_all(intensity_signal)

if measure_intensity == "quantitative: positive nuclei (rec for fluorescence)":
	intensity_signal = quant_intens_positive(imp_signal, positive_nuclear_area, count_nuclei, overlay)	
	show_intensity_tab_positive(intensity_signal)
	

	