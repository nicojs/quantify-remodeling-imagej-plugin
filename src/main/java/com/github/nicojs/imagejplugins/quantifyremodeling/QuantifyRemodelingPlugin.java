package com.github.nicojs.imagejplugins.quantifyremodeling;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import net.imagej.DatasetService;
import org.scijava.ItemVisibility;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

@Plugin(type = Command.class,
        menuPath = "Karin Jansen>Preview Checkbox")
public class QuantifyRemodelingPlugin implements Command, Previewable {

    // -- Parameters --

    @Parameter
    private DatasetService datasetService = null;

    @Parameter(
            autoFill = false,
            description = "Original file description",
            label = "Original file"
    )
    private ImagePlus originalFile;

    @Parameter(
            autoFill = false,
            description = "Orientation file description",
            label = "Orientation file"
    )
    private ImagePlus orientationFile;


    @Parameter(
            autoFill = false,
            description = "Coherency file description",
            label = "Coherency file"
    )
    private ImagePlus coherencyFile;

    @Parameter(
            label = "Border (in microns)"
    )
    private double border;

    @Parameter(
            label = "Offset (in microns)"
    )
    private double offset;

    @Parameter(
            label = "Coherency cutoff (value between 0 and 1)"
    )
    private double coherencyCutoff;

    // -- Command methods --

    @Override
    public void run() {
        double originalWidthPixels = orientationFile.getWidth();
        double originalHeightPixels = originalFile.getHeight();


        if (RoiManager.getInstance() == null) {
            new RoiManager();
        }
        RoiManager roiManager = RoiManager.getInstance();
        roiManager.reset();
        // Using `new Roi()` is equivalent to calling `Roi()` in Jython
        List<Roi> rois = Arrays.asList(
                new Roi(toPixels(offset), 0, toPixels(border), originalHeightPixels),
                new Roi(toPixels(offset + border), 0, toPixels(border), originalHeightPixels),
                new Roi(toPixels(offset + border * 2), 0, toPixels(border), originalHeightPixels)
        );
        rois.forEach(roiManager::addRoi);

        final List<OrderParameter> orderParameters = rois.stream().map(roi -> {
            originalFile.setRoi(roi, true);
            double meansIntensity = originalFile.getStatistics(Measurements.MEAN).mean;
            coherencyFile.setRoi(roi);
            orientationFile.setRoi(roi);


            ImagePlus coherencyRegionImage = coherencyFile.crop();
            ImagePlus orientationRegionImage = orientationFile.crop();
            OrderParameter orderParameter = getOrderParameter(orientationRegionImage, coherencyRegionImage, meansIntensity);
            coherencyRegionImage.close();
            orientationRegionImage.close();
            return orderParameter;
        }).collect(Collectors.toList());

        // https://imagej.net/Tips_for_developers#Show_a_results_table
        ResultsTable orderParameterResults = new ResultsTable();
        for (int i = 0; i < orderParameters.size(); i++) {
            OrderParameter orderParameter = orderParameters.get(i);
            orderParameterResults.incrementCounter(); // next row
            orderParameterResults.addValue("pixel", i);
            orderParameterResults.addValue("um", toMicrons(i));
            orderParameterResults.addValue("param", orderParameter.order_param);
            orderParameterResults.addValue("param_weighted", orderParameter.weighted_order_param);
        }
        orderParameterResults.show("order_param_results.csv");
//
//        for (int i = 0; i < orderParameters.size(); i++) {
//            showBorderResults(orderParameters.get(i), "angle_distribution_" + i + "_borders.csv");
//        }
    }

    private void showBorderResults(OrderParameter orderParameter, String windowTitle) {
        ResultsTable rt = new ResultsTable();
        for (int j = 0; j < orderParameter.coh_list.size(); j++) {
            rt.incrementCounter();
            rt.addValue(orderParameter.or_list.get(j).toString(), orderParameter.coh_list.get(j));
        }
        rt.show(windowTitle);
    }

    private OrderParameter getOrderParameter(ImagePlus orientationRegionImage, ImagePlus coherencyRegionImage, double meansIntensity) {
        float[] orientArray = (float[]) orientationRegionImage.getProcessor().getPixels();
        float[] cohArray = (float[]) coherencyRegionImage.getProcessor().getPixels();

        DoubleStream coherencyStream = IntStream.range(0, cohArray.length)
                .mapToDouble(i -> cohArray[i]);

        List<Double> sin_or = new ArrayList<>();
        List<Double> cos_or = new ArrayList<>();
        List<Double> weight_sin_or = new ArrayList<>();
        List<Double> weight_cos_or = new ArrayList<>();
        List<Double> coh_list = new ArrayList<>();
        List<Double> or_list = new ArrayList<>();

        List<Boolean> mask = coherencyStream
                .mapToObj(x -> x >= coherencyCutoff)
                .collect(Collectors.toList());
        for (int i = 0; i < mask.size(); i++) {
            double sin_tmp = Math.sin(orientArray[i] * 2); // here I do sin(2*angle) where angle is in radian
            double cos_tmp = Math.cos(orientArray[i] * 2);
            double weight_sin_or_tmp = cohArray[i] * Math.sin(orientArray[i] * 2);
            double weight_cos_or_tmp = cohArray[i] * Math.cos(orientArray[i] * 2);

            cos_or.add(cos_tmp);
            sin_or.add(sin_tmp);
            weight_sin_or.add(weight_sin_or_tmp);
            weight_cos_or.add(weight_cos_or_tmp);
            coh_list.add((double) cohArray[i]);
            or_list.add(orientArray[i] * 360 / (2 * Math.PI));
        }

        double avg_sin_or = sin_or.stream().collect(Collectors.averagingDouble(i -> i));
        double avg_cos_or = cos_or.stream().collect(Collectors.averagingDouble(i -> i));

        double avg_weight_sin_or = weight_sin_or.stream().mapToDouble(i -> i).sum() / coh_list.stream().mapToDouble(i -> i).sum();
        double avg_weight_cos_or = weight_cos_or.stream().mapToDouble(i -> i).sum() / coh_list.stream().mapToDouble(i -> i).sum();

        double order_param = Math.sqrt(Math.pow(avg_sin_or, 2) + Math.pow(avg_cos_or, 2));
        double weighted_order_param = Math.sqrt(Math.pow(avg_weight_sin_or, 2) + Math.pow(avg_weight_cos_or, 2));
        return new OrderParameter(order_param, weighted_order_param, coh_list, or_list, meansIntensity);
    }

    private static class OrderParameter {

        private final double meansIntensity;
        private final double order_param;
        private final double weighted_order_param;
        private final List<Double> coh_list;
        private final List<Double> or_list;

        private OrderParameter(double order_param, double weighted_order_param, List<Double> coh_list, List<Double> or_list, double meansIntensity) {
            this.order_param = order_param;
            this.weighted_order_param = weighted_order_param;
            this.coh_list = coh_list;
            this.or_list = or_list;
            this.meansIntensity = meansIntensity;
        }

        private double averageOrientation() {
            return or_list.stream().mapToDouble(i -> i).average().getAsDouble();
        }

        private String toCvsString() {
            NumberFormat englishNumberFormat = NumberFormat.getInstance(Locale.ENGLISH);
            return String.format("%s,%s,%s,%s,%s", "Roi X",
                    englishNumberFormat.format(order_param),
                    englishNumberFormat.format(weighted_order_param),
                    englishNumberFormat.format(meansIntensity),
                    englishNumberFormat.format(averageOrientation()));
        }
    }

    /*

	sin_or=[]
	cos_or=[]
	weight_sin_or=[]
	weight_cos_or=[]
	coh_list=[]
	or_list=[]

	for i in range(len(mask)):
		if mask[i]==True:
			sin_tmp= math.sin(float(orient_array[i])*2) # here I do sin(2*angle) where angle is in radian
			cos_tmp=math.cos(float(orient_array[i]) *2)
			weight_sin_or_tmp=float(coher_array[i])*math.sin(float(orient_array[i])*2)
			weight_cos_or_tmp=float(coher_array[i])*math.cos(float(orient_array[i])*2)
			cos_or.append(cos_tmp)
			sin_or.append(sin_tmp)
			weight_sin_or.append(weight_sin_or_tmp)
			weight_cos_or.append(weight_cos_or_tmp)
			coh_list.append(coher_array[i])
			or_list.append(orient_array[i]*360/(2*math.pi))

	avg_sin_or=float(sum(sin_or)) / float(len(sin_or))
	avg_cos_or=float(sum(cos_or))/float(len(cos_or))
	avg_weight_sin_or=float( sum(weight_sin_or))/ float(sum(coh_list))
	avg_weight_cos_or=float( sum(weight_cos_or))/ float(sum(coh_list))

	order_param=math.sqrt(avg_sin_or**2+avg_cos_or**2)
	weighted_order_param=math.sqrt(avg_weight_sin_or**2+avg_weight_cos_or**2)
	#print (order_param,weighted_order_param)

	return(order_param,weighted_order_param,coh_list,or_list)

    * */


    private double toPixels(double microns) {
        return 1 / originalFile.getCalibration().pixelWidth * microns;
    }

    private double toMicrons(int pixel) {
        return originalFile.getCalibration().pixelWidth * pixel;
    }


    // -- Previewable methods --

    @Override
    public void preview() {

    }

    @Override
    public void cancel() {
        // Set the image's title back to the original value.

    }

    // -- Initializer methods --
    public boolean validateBorder(int value) {
        double maxWidth = orientationFile.getCalibration().pixelWidth * originalFile.getWidth();
        if (value < 0) {
            return false;
        } else return !(value > maxWidth);
    }


}
