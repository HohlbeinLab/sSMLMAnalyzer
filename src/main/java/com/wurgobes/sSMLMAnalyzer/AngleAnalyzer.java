package com.wurgobes.sSMLMAnalyzer;

/*
Spectral Super Resolution Pair Finder
(c) 2021 Martijn Gobes, Wageningen University.

This program takes a list of localisation and using FFT finds features (pairs)
It determines the distance(s) between the pairs and the angle at which they might be positioned

This is far from an infallible process, which is why the user can override this.
Sometimes the result might be flipped n*90 degrees, or mirrored, which is why this script might be run up to 4 times.

This entire process can be replicated in ImageJ to allow for manual tweaks
This is more accurate, but more complicated ofcourse.

The steps to do so are included in the manualAngleAnalyzer.md file

This software is released under the GPL v3. You may copy, distribute and modify
the software as long as you track changes/dates in source files. Any
modifications to or software including (via compiler) GPL-licensed code
must also be made available under the GPL along with build & install instructions.
https://www.gnu.org/licenses/gpl-3.0.en.html

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

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

import java.util.ArrayList;

import static com.wurgobes.sSMLMAnalyzer.Util.*;

import static ij.plugin.filter.ParticleAnalyzer.*;
import static net.imglib2.algorithm.stats.ComputeMinMax.computeMinMax;

public class AngleAnalyzer < T extends IntegerType<T>> implements Command {

    private final LogService logService; // Passed through by instantiated class

    private final FloatMatrix data; // All localisations

    // Results of angles and ditances
    private double angle_low;
    private double angle_high;

    private double dist_low;
    private double dist_high;

    // The settings as set
    private final boolean flipAngles;
    private final boolean mirrorAngles;

    private final boolean debug;

    // Indicates if features were found
    // When there are no regular features (like distances) this will be false
    private boolean succes = false;

    public AngleAnalyzer(FloatMatrix data, boolean flipAngles, boolean mirrorAngles, LogService logService, boolean debug){
        this.data = data;
        this.flipAngles = flipAngles;
        this.mirrorAngles = mirrorAngles;
        this.logService = logService;
        this.debug = debug;

    }

    public void run() {
        logService.info("Analysing Angles and Distances");

        double processingTime = System.nanoTime();

        // FFT works best with images of 2^n size
        // I choose 1024 because it gives good results while still giving good resolution
        // Precision is not a huge issue since we want a ballpark
        // The value is refined later on anyhow
        int target_size = 1024;
        int[] dimensions = new int[] {target_size, target_size};

        // these are the factors why which to reduce each axis to fit inside the image
        final float[] reduce = new float[] {(float) Math.ceil((data.getColumn(1).max() + 1)/target_size), (float) Math.ceil((data.getColumn(2).max() + 1)/target_size)};

        // We can only use FFT on an image, so we create a 1024x1024 image from the points
        ImagePlus temp = getImageFromPoints(data.getColumns(new int[]{1, 2, 3}), reduce, dimensions[0], dimensions[1]);

        if(debug) temp.show();


        ImagePlus FFTImage = FFT.forward(temp); // First FFT with default settings

        Directionality_ direction = new Directionality_(); // Class to calculate directionality (mirror of default ImageJ functionality)

        if(debug) FFTImage.duplicate().show();

        // Set directionality settings
        // We want to detect in a range from 0 to 180 degrees
        direction.setImagePlus(FFTImage);
        direction.setBinRange(0, 180);


        direction.computeHistograms(); // Calculate results

        double[] results = direction.getFitAnalysis().get(0); // Gives the main detected angle, as well as its statistics

        // Get angle and standard deviation
        double angle =  results[0] - Math.PI/2; // rad
        double std = results[1];

        // Get the Processor and rotate it
        // If we calculated the angle wrong we will be able to detect this a bit later
        ImageProcessor FFTProcessor = FFTImage.getProcessor();
        FFTProcessor.rotate(Math.toDegrees(angle)); //degrees


        // The angle can be rotated 180 degrees
        // we try to guess which side it is by splitting the image in twain
        // The 'starting' side should have a higher intensity, and thus a higher summed signal
        // We crop the image around the center and compare the two halves
        // If we need to flip it, we subtract Pi
        Img<T> image = ImageJFunctions.wrapNumeric(FFTImage);

        FinalInterval leftInterval = new FinalInterval(new long[] {256, 256}, new long[] {511, 767});
        FinalInterval rightInterval = new FinalInterval(new long[] {512, 256}, new long[] {767, 767});

        IntervalView<T> left = Views.interval(image, leftInterval);
        IntervalView<T> right = Views.interval(image, rightInterval);

        double different = getSum(left) - getSum(right);

        angle = different < 0 ? angle  : angle - Math.PI;

        if(debug) System.out.println(angle);
        if(debug) System.out.println(std);

        // We rotate the angle by 90 deg
        // or flip it around the x-axis if need be
        if(flipAngles) {
            if (angle > 0) angle -= Math.PI;
            if (angle < 0) angle += Math.PI;
        }
        if(mirrorAngles){
            angle *= -1;
        }

        // The std cant realistically be too high, since they are all split using the same grating
        // if it is too high, we could detect nonsense
        // if it is too low, we could not detect enough, so we set a max and min of 0.2 and 0.04 (rad)
        // Experimentally i have had results have a range of ~0.2 (~0.5 std)
        if(std > 0.2){
            logService.info("Standard Deviation (" + std + ") seems really high." +
                    " Lowering to 0.2 ");
            std = 0.2;
        } else if (std < 0.04) std = 0.04;

        // Set the higher and lower end
        angle_high = angle + 2.5 * std;
        angle_low = angle - 2.5 * std;

        // If either value exceeds Pi, we need to rotate it properly to the other side
        // The angle from atan2() is also calculated in the range [-Pi, Pi)
        if(angle_low < -Math.PI){
            angle_low += 2 * Math.PI;
        }
        if(angle_high > Math.PI){
            angle_high -= 2 * Math.PI;
        }

        // To determine the distance between the pairs we use the feature detection from ImageJ
        // We take only the centre and do another FFT over that
        // The sacrified resolution is not an issue
        FinalInterval centerInterval = new FinalInterval(new long[] {256, 256}, new long[] {767, 767});
        IntervalView<T> center = Views.interval(image, centerInterval);

        ImagePlus centerImgP = ImageJFunctions.wrap(center, "center");
        ImagePlus distance = FFT.forward(centerImgP); // Do FFT

        if(debug) distance.show();

        Img<T> distanceHalf = ImageJFunctions.wrapNumeric(distance);

        // The distance is mirrored around the centre, so we take a small slice from the left to the centre
        // If everything went right the features should show up in this area
        FinalInterval leftDistance = new FinalInterval(new long[] {0, 232}, new long[] {256 ,281});
        IntervalView<T> distanceCenter = Views.interval(distanceHalf, leftDistance);

        ImageJFunctions.show(distanceCenter);

        // We want to filter the image to only the top 0.014% of highest valued pixels
        // To do this we take calculate the highest and lowest value, count the pixels and determine the breakpoint
        // to get to that value
        T min = distanceCenter.firstElement().createVariable(); // Get right variable type
        T max = distanceCenter.firstElement().createVariable();

        computeMinMax(distanceCenter, min, max); // get Min and Max

        // Create a histogram out of the cropped image
        Histogram1d<T> hist = new Histogram1d<>(distanceCenter, new Integer1dBinMapper<>(0, 256, false));

        // Give us a min and max to mask the image with
        int minThreshold = getThresholdBin(0.99_86f, hist.toLongArray());
        int maxThreshold = 256;

        ImagePlus masked = ImageJFunctions.wrap(distanceCenter, "masked"); // create a copy of the image to mask

        if(debug)  masked.show();

        masked.getProcessor().setThreshold(minThreshold, maxThreshold, ImageProcessor.BLACK_AND_WHITE_LUT); // set threshhold

        (new Thresholder()).run("mask"); // execute threshold

        ResultsTable resultsTable = new ResultsTable(); // get a resultstable to put the feature finding results into

        // Setup the feature analyzer to get the center of mass of any blobs
        ParticleAnalyzer particleAnalyzer = new ParticleAnalyzer(SHOW_MASKS, Analyzer.getMeasurements() | Analyzer.CENTER_OF_MASS, resultsTable, 0, 100);
        particleAnalyzer.analyze(masked);

        if(!debug) WindowManager.closeAllWindows();

        // get the results from the table
        float[] offsetsX = resultsTable.getColumn(resultsTable.getColumnIndex("XM"));
        float[] offsetsY = resultsTable.getColumn(resultsTable.getColumnIndex("YM"));

        // sort them so the first value is the centre (and thus most right) and the next values are left from that
        float[][] sortResult = sortMultiple(offsetsX, offsetsY);


        offsetsX = sortResult[0];
        offsetsY = sortResult[1];

        // each feature (minus the centre) is a peak
        int peaks = offsetsX.length - 1;

        // This is the scale of units/px
        // both the original and FFT are target size
        float sizeX = data.getColumn(1).max() / 512; // deltaW
        float sizeY = data.getColumn(2).max() / 512;


        // Analyze the results, discarding any that are too far away or too close
        // this also calculates the actual distance to the centre and its angle
        // if the angle is > 0.1 we did not do a good calculation and echo this back
        // as well as correct it using this check
        ArrayList<Double> angles = new ArrayList<>();
        ArrayList<Double> distances = new ArrayList<>();

        // These cutoff values are chosen empirically
        // With the FFT these distances can never vary a huge amount
        // If it turns out someone finds an application where a larger size can occur, one can change these
        float averageSize = (data.getColumn(1).max() + data.getColumn(2).max())/2;
        double lowerCutoff = 0.01 * averageSize; // Distances between features is at min 1% of whole size
        double upperCutoff = 0.1 * averageSize; // Distances between features is at most 10% of whole size

        for(int i = 1; i <= peaks; i++){
            float posX = offsetsX[0] - offsetsX[i]; // delta x to centre
            float posY = offsetsY[0] - offsetsY[i]; // delta y to centre

            // calculate distance to center of feature
            double dist = Math.sqrt(Math.pow(posX * sizeX, 2) + Math.pow(posY * sizeY, 2));

            // If the feature is too close or too far, it probably is an error (and would give really wrong results
            if(dist < lowerCutoff | dist > upperCutoff) continue;

            distances.add(dist);

            double curAngle = Math.atan2(posY, posX);
            if(debug) System.out.println(curAngle);
            if(Math.abs(curAngle) > 0.1){ // Check angle, shouldnt be far from 0
                angles.add(curAngle);
            }
        }

        // Correct the angle if we still detected one in the feature finding
        if(angles.size() > 1) {
            double realAngle = angles.stream()
                    .mapToDouble(a -> a)
                    .average().orElse(Double.NaN);

            angle_low = angle - realAngle - std*2;
            angle_high = angle - realAngle + std*2;
            if(debug)
                logService.info("Calculating the angle was inaccurate. Angle seems to be: " + realAngle);
        }

        // If we find features we had succes!
        // The actual distance is calculated with some extra buffer to ensure we get all points
        // we have no std due to low sample size so these have a mostly emperical adjustment
        // the low end is given by the first peak
        // the high end by the last peak
        if(distances.size() > 0){
            succes = true;
            double buffer = 0.0075 * averageSize;

            dist_low = distances.get(0) * 0.90 - buffer/2;
            dist_high = distances.get(distances.size() - 1)*1.10 + buffer/2;

            logService.info("Angle: " + angle_low  + "-" + angle_high );
            logService.info("Distance: " + Math.round(dist_low) + "-" + Math.round(dist_high));
        }


        processingTime = System.nanoTime() - processingTime;
        System.out.println("Calculating Angles and Distances took " + String.format("%.3f", processingTime / 1000000000) + " s");
    }

    public float[] getAngles(){
        //
        return new float[]{(float) angle_low, (float)  angle_high };
    }

    public float[] getDistances(){
        return new float[]{Math.round(dist_low), Math.round(dist_high)};
    }

    public boolean getSucces() {return succes;}

}
