package com.wurgobes.sSMLMAnalyzer;


import fiji.analyze.directionality.Directionality_;
import ij.*;

import ij.measure.ResultsTable;
import ij.plugin.FFT;

import ij.plugin.Thresholder;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;

import net.imglib2.FinalInterval;


import net.imglib2.histogram.Histogram1d;
import net.imglib2.histogram.Integer1dBinMapper;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.IntegerType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import org.jblas.FloatMatrix;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Plugin;



import java.util.ArrayList;

import static com.wurgobes.sSMLMAnalyzer.Util.*;

import static ij.plugin.filter.ParticleAnalyzer.*;
import static net.imglib2.algorithm.stats.ComputeMinMax.computeMinMax;

@Plugin(type = Command.class, menuPath = "Plugins>Spectral Analyzer>Analyze Angles and Distance")
public class AngleAnalyzer < T extends IntegerType<T>> implements Command {

    private final LogService logService;


    private final FloatMatrix data;

    private double angle_low;
    private double angle_high;

    private double dist_low;
    private double dist_high;

    private boolean flipAngles = false;

    public AngleAnalyzer(FloatMatrix data, boolean flipAngles, LogService logService){
        this.data = data;
        this.flipAngles = flipAngles;
        this.logService = logService;
    }

    @Override
    public void run() {
        logService.info("Analysing Angles and Distances");

        double processingTime = System.nanoTime();

        // Implement FFT to guess angle and distance

        int target_size = 1024;
        int[] dimensions = new int[] {target_size, target_size};
        final float[] reduce = new float[] {(float) Math.ceil((data.getColumn(1).max() + 1)/target_size), (float) Math.ceil((data.getColumn(2).max() + 1)/target_size)};

        ImagePlus temp = getImageFromPoints(data.getColumns(new int[]{1, 2, 3}), reduce, dimensions[0], dimensions[1]);

        //temp.show();

        ImagePlus FFTImage = FFT.forward(temp);

        Directionality_ direction = new Directionality_();

        //FFTImage.duplicate().show();

        direction.setImagePlus(FFTImage);
        direction.setBinRange(0, 180);

        direction.computeHistograms();

        double[] results = direction.getFitAnalysis().get(0);

        double angle =  results[0] - Math.PI/2;
        double std = results[1];


        ImageProcessor FFTProcessor = FFTImage.getProcessor();
        FFTProcessor.rotate(Math.toDegrees(angle)); //degrees


        Img<T> image = ImageJFunctions.wrapNumeric(FFTImage);


        FinalInterval leftInterval = new FinalInterval(new long[] {256, 256}, new long[] {511, 767});
        FinalInterval rightInterval = new FinalInterval(new long[] {512, 256}, new long[] {767, 767});

        IntervalView<T> left = Views.interval(image, leftInterval);
        IntervalView<T> right = Views.interval(image, rightInterval);

        double different = getSum(left) - getSum(right);


        angle = different < 0 ? angle  : angle - Math.PI;

        if(flipAngles) angle -= Math.PI;

        if(std > 0.2){
            logService.info("Standard Deviation (" + std + ") seems really high." +
                    " Lowering to 0.2 ");
            std = 0.2;
        } else if (std < 0.04) std = 0.04;

        angle_high = angle + 1.5 * std;
        angle_low = angle - 1.5 * std;

        if(angle_low < -Math.PI){
            angle_low += 2 * Math.PI;
        }
        if(angle_high > Math.PI){
            angle_high -= 2 * Math.PI;
        }

        FinalInterval centerInterval = new FinalInterval(new long[] {256, 256}, new long[] {767, 767});
        IntervalView<T> center = Views.interval(image, centerInterval);

        ImagePlus centerImgP = ImageJFunctions.wrap(center, "center");
        ImagePlus distance = FFT.forward(centerImgP);

        //distance.show();

        Img<T> distanceHalf = ImageJFunctions.wrapNumeric(distance);

        FinalInterval leftDistance = new FinalInterval(new long[] {0, 232}, new long[] {256 ,281});
        IntervalView<T> distanceCenter = Views.interval(distanceHalf, leftDistance);

        ImageJFunctions.show(distanceCenter);

        T min = distanceCenter.firstElement().createVariable();
        T max = distanceCenter.firstElement().createVariable();

        computeMinMax( distanceCenter, min, max );

        Histogram1d<T> hist = new Histogram1d<>(distanceCenter, new Integer1dBinMapper<>(0, 256, false));

        int minThreshold = getThresholdBin(0.99_86f, hist.toLongArray());
        int maxThreshold = 256;

        ImagePlus masked = ImageJFunctions.wrap(distanceCenter, "masked");

        //masked.show();

        masked.getProcessor().setThreshold(minThreshold, maxThreshold, ImageProcessor.BLACK_AND_WHITE_LUT);

        (new Thresholder()).run("mask");

        ResultsTable resultsTable = new ResultsTable();


        ParticleAnalyzer particleAnalyzer = new ParticleAnalyzer(SHOW_MASKS, Analyzer.getMeasurements() | Analyzer.CENTER_OF_MASS, resultsTable, 0, 100);
        particleAnalyzer.analyze(masked);

        WindowManager.closeAllWindows();

        float[] offsetsX = resultsTable.getColumn(resultsTable.getColumnIndex("XM"));
        float[] offsetsY = resultsTable.getColumn(resultsTable.getColumnIndex("YM"));


        float[][] sortResult = sortMultiple(offsetsX, offsetsY);


        offsetsX = sortResult[0];
        offsetsY = sortResult[1];


        int peaks = offsetsX.length - 1;

        // both the original and FFT are target size
        float sizeX = data.getColumn(1).max() / target_size; // deltaW
        float sizeY = data.getColumn(2).max() / target_size;

        double constMul = 2.09; // Got this emperically, want to find a better value (with logic behind it)



        ArrayList<Double> angles = new ArrayList<>();
        ArrayList<Double> distances = new ArrayList<>();
        for(int i = 1; i <= peaks; i++){
            float posX = offsetsX[0] - offsetsX[i];
            float posY = offsetsY[0] - offsetsY[i];

            double dist = constMul * Math.sqrt(Math.pow(posX * sizeX, 2) + Math.pow(posY * sizeY, 2));

            if(dist < 500 | dist > 4000)
                continue;


            distances.add(dist);

            double curAngle = Math.atan2(posY, posX);
            if(Math.abs(curAngle) > 0.1){
                angles.add(curAngle);
            }
        }

        if(angles.size() > 1) {
            double realAngle = angles.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(Double.NaN);

            angle_low = realAngle - std*0.75;
            angle_high = realAngle + std*0.75;
            logService.info("Calculating the angle was inaccurate. Angle seems to be: "
                    + realAngle
            );
        }



        double buffer = 250;

        dist_low = distances.get(0) - buffer/2;
        dist_high = distances.get(distances.size() - 1) + buffer/2;

        logService.info("Angle: " + angle_low  + "-" + angle_high );
        logService.info("Distance: " + Math.round(dist_low) + "-" + Math.round(dist_high));


        processingTime = System.nanoTime() - processingTime;
        System.out.println("Calculating Angles and Distances took " + String.format("%.3f", processingTime / 1000000000) + " s");
    }

    public float[] getAngles(){
        return new float[]{(float)  angle_low, (float)  angle_high };
    }

    public float[] getDistances(){
        return new float[]{(float) Math.round(dist_low), (float) Math.round(dist_high)};
    }

    public static void main(String[] args) {
        //Change data to be static before using this
        /*
        ImageJ ij = new ImageJ();
        ij.ui().showUI();


        final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};


        String filePath = "F:\\ThesisData\\output\\4_grating_drift.csv";
        // String filePath = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\temp.csv";
        // String filePath = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\distance_testing.csv";

        FloatMatrix floatMatrix = null;
        List<String> collumns = new ArrayList<>();

        OwnFloatMatrixLoader ownFloatMatrix = new OwnFloatMatrixLoader();

        try {
            floatMatrix =  ownFloatMatrix.loadCSVFile(filePath);
            collumns = ownFloatMatrix.collumns;
        } catch (IOException e) {
            System.out.println("File not found.");
        } catch (LapackException e) {
            e.printStackTrace();
        } catch (Exception e){
            System.out.println("Does the csv start with a header?");
            e.printStackTrace();
        }

        assert !collumns.isEmpty();
        assert floatMatrix != null;


        int[] revOptionsIndices = new int[possible_options.length];


        Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

        for(int i = 0; i < collumns.size(); i ++){
            String header = collumns.get(i);
            Matcher matcher = pattern.matcher(header);
            if(matcher.find()){
                revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;
            }
        }

        //frame, x, y, intensity
        data2 = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3], revOptionsIndices[5]});

        ij.command().run(AngleAnalyzer.class, true);

         */

    }
}
