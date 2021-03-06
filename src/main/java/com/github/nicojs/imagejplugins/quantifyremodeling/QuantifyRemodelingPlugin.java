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
        menuPath = "Plugins>Quantify remodeling")
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
            label = "Border (in microns)",
            validater = "validateBorder"
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
            orderParameterResults.addValue("Window nr.", i + 1);
            orderParameterResults.addValue("Window size (um)", border);
            orderParameterResults.addValue("param", orderParameter.order_param);
            orderParameterResults.addValue("param_weighted", orderParameter.weighted_order_param);
            orderParameterResults.addValue("fluorescent_Intensity (a.u.)", orderParameter.meansIntensity);
            orderParameterResults.addValue("average_orientation", orderParameter.average_orientation);
        }
        orderParameterResults.show("Results Intensity and order parameters");

        for (int i = 0; i < orderParameters.size(); i++) {
            showBorderResults(orderParameters.get(i), "ROI_" + (i + 1) + "_angle_distribution");
        }
    }

    private void showBorderResults(OrderParameter orderParameter, String windowTitle) {
        ResultsTable rt = new ResultsTable();
        for (int j = 0; j < orderParameter.coh_list.size(); j++) {
            rt.incrementCounter();
            rt.addValue("Orientation", orderParameter.or_list.get(j).toString());
            rt.addValue("Coherence", orderParameter.coh_list.get(j));
        }
        rt.show(windowTitle);
    }

    private OrderParameter getOrderParameter(ImagePlus orientationRegionImage, ImagePlus coherencyRegionImage, double meansIntensity) {
        float[] orientationPixels = (float[]) orientationRegionImage.getProcessor().getPixels();
        float[] coherencyPixels = (float[]) coherencyRegionImage.getProcessor().getPixels();

        List<Double> sin_or = new ArrayList<>();
        List<Double> cos_or = new ArrayList<>();
        List<Double> weight_sin_or = new ArrayList<>();
        List<Double> weight_cos_or = new ArrayList<>();
        List<Double> coh_list = new ArrayList<>();
        List<Double> or_list = new ArrayList<>();

        for (int i = 0; i < coherencyPixels.length; i++) {
            double coherencyValue = coherencyPixels[i];
            if (coherencyValue >= coherencyCutoff) {
                double sin_tmp = Math.sin(orientationPixels[i] * 2); // here I do sin(2*angle) where angle is in radian
                double cos_tmp = Math.cos(orientationPixels[i] * 2);
                double weight_sin_or_tmp = coherencyValue * sin_tmp;
                double weight_cos_or_tmp = coherencyValue * cos_tmp;

                cos_or.add(cos_tmp);
                sin_or.add(sin_tmp);
                weight_sin_or.add(weight_sin_or_tmp);
                weight_cos_or.add(weight_cos_or_tmp);
                coh_list.add(coherencyValue);
                or_list.add(toDegree(orientationPixels[i]));
            }
        }

        double avg_sin_or = sin_or.stream().collect(Collectors.averagingDouble(i -> i));
        double avg_cos_or = cos_or.stream().collect(Collectors.averagingDouble(i -> i));

        final double totalCoherency = coh_list.stream().mapToDouble(i -> i).sum();
        // Take the weighted average of sin and cos
        double avg_weight_sin_or = weight_sin_or.stream().mapToDouble(i -> i).sum() / totalCoherency;
        double avg_weight_cos_or = weight_cos_or.stream().mapToDouble(i -> i).sum() / totalCoherency;

        double order_param = Math.sqrt(Math.pow(avg_sin_or, 2) + Math.pow(avg_cos_or, 2));
        double weighted_order_param = Math.sqrt(Math.pow(avg_weight_sin_or, 2) + Math.pow(avg_weight_cos_or, 2));
        double average_orientation = or_list.stream().collect(Collectors.averagingDouble(i -> i));
        return new OrderParameter(order_param, weighted_order_param, coh_list, or_list, meansIntensity, average_orientation);
    }

    private double toDegree(double orientationPixel) {
         return orientationPixel * 360 / (2 * Math.PI);
    }

    private static class OrderParameter {
        //need to add average_orientation in here as well
        private final double meansIntensity;
        private final double order_param;
        private final double weighted_order_param;
        private final List<Double> coh_list;
        private final List<Double> or_list;
        private final double average_orientation;

        private OrderParameter(double order_param, double weighted_order_param, List<Double> coh_list, List<Double> or_list, double meansIntensity, double average_orientation) {
            this.order_param = order_param;
            this.weighted_order_param = weighted_order_param;
            this.coh_list = coh_list;
            this.or_list = or_list;
            this.meansIntensity = meansIntensity;
            this.average_orientation = average_orientation;
        }
    }


    private int toPixels(double microns) {
        return (int) Math.round(1 / originalFile.getCalibration().pixelWidth * microns);
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
    public boolean validateBorder() {
//        double maxWidth = orientationFile.getCalibration().pixelWidth * originalFile.getWidth();
//        if (value < 0) {
//            return false;
//        } else return !(value > maxWidth);
        System.out.println("test");
        return true;
    }

}
