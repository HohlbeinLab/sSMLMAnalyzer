package com.wurgobes.sSMLMAnalyzer;

/*
Spectral Super Resolution Pair Finder
(c) 2021 Martijn Gobes, Wageningen University.

This file contains various utility functions used by the other files

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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import net.imglib2.type.numeric.RealType;
import net.imglib2.util.RealSum;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jblas.FloatMatrix;
import org.jblas.ranges.IntervalRange;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import static ij.util.ThreadUtil.createThreadArray;
import static ij.util.ThreadUtil.startAndJoin;
import static org.jblas.MatrixFunctions.*;


public class Util {

    public static FloatMatrix atan2(FloatMatrix x, FloatMatrix y){
        // Short helper function for atan2
        return atan2(x, y, Distance(x, y));
    }

    public static FloatMatrix atan2(final FloatMatrix x, final FloatMatrix y, final FloatMatrix intermediate){
        // Calculated the four quadrant inverse tangent for the matrix with x and y
        //Very much not a perfect solution
        // Method taken from https://en.wikipedia.org/wiki/Atan2
        x.assertSameSize(y);
        final FloatMatrix result = new FloatMatrix(x.rows, x.columns);

        final FloatMatrix domainOne = atan(y.div(intermediate.add(x))).mul(2);

        for(int index : x.gt(0).or(y.ne(0)).findIndices()) result.put(index, domainOne.get(index));
        for(int index : x.lt(0).and(y.eq(0)).findIndices()) result.put(index, (float) Math.PI);
        for(int i = 0; i < result.rows; i++) result.put(i, i, Float.NaN);

        return result;
    }

    public static FloatMatrix Distance(final FloatMatrix X, final FloatMatrix Y){
        // Calculates distance from each point to each other point
        return sqrt(X.mul(X).addi(Y.mul(Y)));
    }

    public static FloatMatrix makeSubstractedMatrix(FloatMatrix A){
        // Creates matrix that shows the distance from one value to each other value (including itself)

        if(A.columns > 1)
            throw new IllegalArgumentException("Matrix can have only 1 collumn");

        int dataPoints = A.rows;

        FloatMatrix subtracted = new FloatMatrix(dataPoints, dataPoints);

        for(int point = 0; point < dataPoints; point++){
            float a = A.get(point);
            subtracted.putRow(point, A.rsub(a));
        }
        return subtracted;
    }

    // Some functions that create a double array because they're needed for ImageJ functions
    public static double[] toDouble(FloatMatrix A){
        double[] result = new double[A.length];

        for(int i = 0; i < A.length; i++) result[i] = A.get(i);

        return result;
    }

    public static double[] toDouble(float v){
        return new double[] {v};
    }

    public static double[] toDouble(List<Integer> list){
        return list.stream().mapToDouble(i->i).toArray();
    }

    public static FloatMatrix extend(final FloatMatrix A, int rows, int columns){
        // extends a FloatMatrix with 0's to the size provided
        if(A.rows == rows && A.columns == columns) return A;
        FloatMatrix result = new FloatMatrix(rows, columns);
        for(int i = 0; i < A.rows; i++){
            for(int j = 0; j < A.columns; j++){
                result.put(i, j, A.get(i, j));
            }
        }
        return result;
    }

    public static void SaveCSV(final FloatMatrix data, List<String> Headers, Path CSV_FILE_NAME)  {
        // Creates a csv file and writes all the data to it
        File csvOutputFile = CSV_FILE_NAME.toFile();
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            pw.println(String.join(",", Headers));
            pw.println(data.toString("%f", "", "",", ", "\n"));
        } catch(IOException error) {
            System.out.println("Could not save CSV.");
            error.printStackTrace();
        }
    }

    public static int getBins(FloatMatrix A, float width){
        // Returns the amount of bins required to get the width desired
        return (int) ((A.max()-A.min())/width);
    }

    public static String getOrdinal(int i) {
        // Gives you the ordinal associated with the number provided
        String[] suffixes = new String[] { "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th" };
        switch (i % 100) {
            case 11:
            case 12:
            case 13:
                return i + "th";
            default:
                return i + suffixes[i % 10];
        }
    }

    public static String getTitleHist(int i){
        // Creates a nicely formatted title for a histogram
        return(getOrdinal(i) + "-" + getOrdinal(i + 1));
    }

    public static String getTitlePlot(int orders){
        // Creates a nicely formatted title for a plot
        StringBuilder title = new StringBuilder();
        for(int i = 0; i < orders; i++){
            title.append(getOrdinal(i)).append("\n");
        }
        return title.toString();
    }

    public static FloatMatrix cleanup(final FloatMatrix A, int neighbours, float distance, int coreCount){
        // Clean up the matrix by discarding any points that do not have at least N neighbours within D distance of them

        final List<Integer> test = Collections.synchronizedList(new ArrayList<>());

        final FloatMatrix X = A.getColumn(3);
        final FloatMatrix Y = A.getColumn(4);

        final AtomicInteger ai = new AtomicInteger(0);
        final Thread[] threads = createThreadArray(coreCount); //Get maximum of threads

        final int rows = A.rows;
        //Set the run function for each thread
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread(() -> {
                for (int row = ai.getAndIncrement(); row <= rows; row = ai.getAndIncrement()) {
                    final int currentAI = ai.get();
                    if (currentAI % 1000 == 0) System.out.println("\r" + currentAI + "/" + rows);
                    IJ.showProgress(currentAI, rows);
                    IJ.showStatus(currentAI + "/" + rows);

                    final float x = A.get(row, 3);
                    final float y = A.get(row, 4);

                    final FloatMatrix distances = Distance(X.sub(x), Y.sub(y));

                    if(distances.lt(distance).sum() > neighbours + 1) test.add(row);
                }
            });
        }

        startAndJoin(threads);

        final int[] indices =  test.stream().mapToInt(i->i).toArray();
        return A.getRows(indices);
    }

    public static FloatMatrix connectOrders(FloatMatrix intermediate, int orders, int orderColumns){
        // This function takes all points and tried to connect the different points together for at most N orders.
        // If more points are found this is echo'd
        boolean checkMoreOrders = true;
        int max_orders = 5;
        if(orders > max_orders) max_orders = orders;
        final FloatMatrix intermediateExtended = extend(intermediate, intermediate.rows, max_orders * orderColumns);
        List<Integer> toKeep = new ArrayList<>();
        for(int row = 0; row < intermediateExtended.rows; row ++){
            float index = intermediateExtended.get(row, 2);

            if(intermediateExtended.getColumn(7).eq(index).sum() == 0.0f) { //Start of chain
                toKeep.add(row);

                int connected_to = (int) intermediateExtended.get(row, 7);
                int[] connected_indices = intermediateExtended.getColumn(2).eq(connected_to).findIndices();

                if(connected_indices.length > 0){
                    //We already did the 0-1 connection

                    intermediateExtended.putRow(row, recursiveSearch(intermediateExtended, intermediateExtended.getRow(row), connected_indices, 2, max_orders, orderColumns));

                    if(checkMoreOrders){
                        connected_to = (int) intermediateExtended.get(row, (orders-1) * orderColumns);

                        if(intermediateExtended.getColumn(2).eq(connected_to).sum() > 1.0f) {
                            System.out.println("There seem to be more orders than " + orders);
                            checkMoreOrders = false;
                        }
                    }
                }
            }
        }

        int[] indices = toKeep.stream().mapToInt(i->i).toArray();

        return intermediateExtended.getColumns(new IntervalRange(0, orders * orderColumns)).getRows(indices);
    }

    private static FloatMatrix recursiveSearch(FloatMatrix intermediate, FloatMatrix rowData, int[] connected_indices, int order, int max_order, int orderCollumns) {
        // Helper function for Connect orders that connects orders by searching for any possible matches and if multiple are found, tries to find the best
        // Ties are broken by which point has the closest angle to the original pair
        // This is recursive to find all next orders

        //i.e. 14: id, 15: x, 16: y, 17: z, 18: intensity, 19: distance, 20: angle
        int[] target_range = new int[]{order * orderCollumns, 1 + (order * orderCollumns), 2 + (order * orderCollumns), 3 + (order * orderCollumns), 4 + (order * orderCollumns), 5 + (order * orderCollumns), 6 + (order * orderCollumns)};

        if (connected_indices.length == 1) {
            FloatMatrix target_row = intermediate.getRow(connected_indices[0]);
            rowData.put(target_range, target_row.getColumns(new int[]{7, 8, 9, 10, 11, 12, 13}));
        } else {
            FloatMatrix intermediateResults = new FloatMatrix(connected_indices.length, rowData.columns);

            for (int i = 0; i < connected_indices.length; i++) {

                FloatMatrix target_row = intermediate.getRow(connected_indices[i]);
                intermediateResults.putRow(i, rowData.dup().put(target_range, target_row.getColumns(new int[]{7, 8, 9, 10, 11, 12, 13})));
                int[] nextConnectedIndices = intermediate.getColumn(2).eq(target_row.get(target_range[0])).findIndices();
                if (nextConnectedIndices.length > 0 && order + 1 < max_order)
                    intermediateResults.putRow(i, recursiveSearch(intermediate, intermediateResults.getRow(i), nextConnectedIndices, order + 1, max_order, orderCollumns));

            }
            intermediateResults = intermediateResults.getRows(FilterbyLength(intermediateResults, max_order));
            if(intermediateResults.rows == 1){
                return intermediateResults;
            } else{
                float originalAngle = rowData.get(6 + ((order - 1) * orderCollumns));

                return intermediateResults.getRow(getClosestIndex(intermediateResults.getColumn(target_range[4]), originalAngle));
            }
        }

        int[] nextConnectedIndices = intermediate.getColumn(2).eq(rowData.get(target_range[0])).findIndices();

        if (nextConnectedIndices.length == 0 | order + 1 == max_order) {
            return rowData;
        } else {
            return recursiveSearch(intermediate, rowData, nextConnectedIndices, order + 1, max_order, orderCollumns);
        }

    }

    public static int[] FilterbyLength(final FloatMatrix A, int orders){
        // Find if there are any points left in a certain order
        for(int i = orders - 1; i > 1; i--){
            FloatMatrix IdCollumn = A.getColumn(i * 5);
            if(IdCollumn.sum() > 0.0f){
                return IdCollumn.ne(0.0f).findIndices();
            }
        }
        return new int[] {0};
    }

    public static int getClosestIndex(final FloatMatrix A, float v){
        // Find the closest index to value v
        final FloatMatrix minimised = abs(A.sub(v));
        return minimised.eq(minimised.min()).findIndices()[0];
    }

    public static FloatMatrix abs(final FloatMatrix A){
        // do inplace abs calculation on the matrix
        for(int i = 0; i < A.length; i++) A.put(i, Math.abs(A.get(i)));
        return A;
    }

    public static ImagePlus getImageFromPoints(final FloatMatrix A, float[] reduction, int width, int height){
        // Creates an imageplus from the Floatmatrix provided for the given height and width
        // It divides each axis by the value provided since it might be too much otherwise
        // These values are chosen in such a way it is never too large
        float[] data = new float[width * height];
        float max_intensity = A.getColumn(2).max();
        FloatMatrix intensities = A.getColumn(2).div(max_intensity);
        for(int i = 0; i < A.rows; i++){
            data[(int)(A.get(i, 0)/reduction[0]) + (int)(A.get(i, 1)/reduction[1]) * height] = intensities.get(i);
        }

        return new ImagePlus("Points", new FloatProcessor(width, height, data, null));
    }

    public static < T extends RealType< T > > double getSum( final Iterable< T > iterable ) {
        // Sums all values in the iterable
        final RealSum sum = new RealSum();

        for ( final T type : iterable )
            sum.add( type.getRealDouble() );

        return sum.getSum();
    }

    public static void addLutLegend(Plot plot, OwnColorTable ct, String label, int width, double start, double end){
        // Add a LUT legend to the provided plot with the provided start, end and title
        for(int i = 0; i < width; i++){
            plot.setColor(ct.getColor(i, 0, width));
            plot.drawNormalizedLine(0.01 + 0.0005 * (i+2),0.93, 0.01 + 0.0005 * (i+2), 0.99);
        }
        plot.setColor("black");
        plot.drawNormalizedLine(0.01,0.925, 0.01, 0.995);
        plot.drawNormalizedLine(0.01 + 0.0005 * (width+2),0.925, 0.01 + 0.0005 * (width+2), 0.995);
        plot.drawNormalizedLine(0.011, 0.99, 0.01 + 0.0005 * (width+2), 0.99);

        double diff = end - start;

        for(double i = start + 25; i < end; i += 25){
            double w = ((i-start)/diff)*width;
            plot.drawNormalizedLine(0.01 + 0.0005 * (w+2), 0.986, 0.01 + 0.0005 * (w+2), 0.99);
        }

        for(double i = start + 50; i < end; i += 50){
            double w = ((i-start)/diff)*width;
            plot.drawNormalizedLine(0.01 + 0.0005 * (w+2), 0.984, 0.01 + 0.0005 * (w+2), 0.99);
        }

        for(double i = start + 100; i < end; i += 100){
            double w = ((i-start)/diff)*width;
            plot.drawNormalizedLine(0.01 + 0.0005 * (w+2), 0.982, 0.01 + 0.0005 * (w+2), 0.99);
        }

        plot.addLabel(0.001, 0.92, String.valueOf((int) start));
        plot.addLabel(0.001 + 0.00025 * (width-5), 0.92, label);
        plot.addLabel(0.001 + 0.0005 * (width+2), 0.92, String.valueOf((int) end));
    }

    public static long sum(long[] rx) {
        // long sum
        long sum = 0L;

        for (long l : rx) {
            sum += l;
        }

        return sum;
    }

    public static int sum(boolean[] rx){
        // count 'true' items
        int sum = 0;
        for(boolean b : rx){
            sum += b ? 1 : 0;
        }
        return sum;
    }

    public static int getThresholdBin(float threshold, final long[] hist){
        // Get the amount of items that are under the threshold in the histogram
        int bins = hist.length;
        final long elementSum = sum(hist);

        threshold *= elementSum;

        long temp = 0L;
        for(int i = 0; i < bins; i++){
            temp += hist[i];
            if(temp > threshold) return i;
        }

        return 0;
    }

    public static float[][] sortMultiple(float[] a, float[] b){
        // Sorts a, and then also sorts b index by index
        float[] tempA = a.clone();
        float[] tempB = new float[b.length];


        Arrays.sort(tempA);
        reverse(tempA);

        for(int i = 0; i < a.length; i++){
            for(int j = 0; j < a.length; j++){
                if(a[j]==tempA[i]){
                    tempB[i] = b[j];

                    break;
                }
            }
        }
        //return new float[][]{tempA, tempB, tempC};
        return new float[][]{tempA, tempB};
    }

    public static void reverse(final float[] input) {
        // reverses the array
        int last = input.length - 1;
        int middle = input.length / 2;
        for (int i = 0; i <= middle; i++) {
            float temp = input[i];
            input[i] = input[last - i];
            input[last - i] = temp;
        }
    }

    public static boolean[][] checkForRetry(FloatMatrix A){
        // Provide some checks on the list to see fi it had a guassian shape-ish, and if the tailbins contain more values than expected
        StandardDeviation std = new StandardDeviation();

        double min = A.min();
        double max = A.max();

        double width = 0.005f;

        double lowerbins = A.le((float) (min + width)).sum();
        double upperbins = A.ge((float) (max - width)).sum();

        List<Integer> bins = new ArrayList<>();
        for(double i = min; i < max; i += width){
            bins.add((int) A.ge((float) (i - width)).and(A.le((float) (i + width))).sum());
        }

        int maxVal = Collections.max(bins);
        int minVal = Collections.min(bins);
        int maxIdx = bins.indexOf(maxVal);
        int minIdx = bins.indexOf(minVal);

        double buffer = 0.2;

        return new boolean[][] {{(std.evaluate(toDouble(bins)) / A.rows) > 0.05, minIdx < bins.size()*buffer || minIdx > bins.size()*(1-buffer), maxIdx > (bins.size() * 0.4) || maxIdx < (bins.size() * 0.6)}, {Math.abs(upperbins - lowerbins) > Math.max(upperbins, lowerbins) * 0.2, lowerbins > upperbins}};
    }

    public static void saveThunderSTORM(Path CSV_FILE_NAME, final FloatMatrix data){
        List<String> ShortHeader = new ArrayList<>();
        ShortHeader.add("id");
        ShortHeader.add("frame");
        ShortHeader.add("x [nm]");
        ShortHeader.add("y [nm]");
        ShortHeader.add("intensity [photons]");
        ShortHeader.add("z [nm]");
        // Creates a csv file and writes all the data to it
        File csvOutputFile = CSV_FILE_NAME.toFile();
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            pw.println(String.join(",", ShortHeader));
            pw.println(data.toString("%f", "", "",", ", "\n"));
        } catch(IOException error) {
            System.out.println("Could not save CSV.");
            error.printStackTrace();
        }
    }

    public static FloatMatrix sort(FloatMatrix matrix, int column){
        // Sorts the matrix based on a column
        FloatMatrix result = new FloatMatrix(matrix.rows, matrix.columns);

        int[] permutation = matrix.getColumn(column).sortingPermutation();

        for(int i = 0; i < matrix.rows; i++){
            result.putRow(i, matrix.getRow(permutation[i]));
        }

        return result;
    }
}
