package com.wurgobes.sSMLMAnalyzer;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import net.imglib2.type.numeric.RealType;
import net.imglib2.util.RealSum;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.jblas.FloatMatrix;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
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
        return atan2(x, y, Distance(x, y));
    }

    public static FloatMatrix atan2(FloatMatrix x, FloatMatrix y, FloatMatrix intermediate){
        //Very much not a perfect solution
        //https://en.wikipedia.org/wiki/Atan2
        x.assertSameSize(y);
        FloatMatrix result = new FloatMatrix(x.rows, x.columns);

        FloatMatrix domainOne = atan(y.div(intermediate.add(x))).mul(2);
        //FloatMatrix domainTwo = atan(intermediate.sub(x).div(y)).mul(2);

        for(int index : x.gt(0).or(y.ne(0)).findIndices()) result.put(index, domainOne.get(index));
        //for(int index : x.le(0).and(y.ne(0)).findIndices()) result.put(index, domainTwo.get(index));
        for(int index : x.lt(0).and(y.eq(0)).findIndices()) result.put(index, (float) Math.PI);
        for(int i = 0; i < result.rows; i++) result.put(i, i, Float.NaN);


        return result;
    }

    public static FloatMatrix Distance(FloatMatrix X, FloatMatrix Y){ return sqrt(X.mul(X).addi(Y.mul(Y))); }

    public static FloatMatrix makeSubstractedMatrix(FloatMatrix A){
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

    public static FloatMatrix extend(FloatMatrix A, int rows, int collumns){
        FloatMatrix result = new FloatMatrix(rows, collumns);
        for(int i = 0; i < A.rows; i++){
            for(int j = 0; j < A.columns; j++){
                result.put(i, j, A.get(i, j));
            }
        }
        return result;
    }

    public static void SaveCSV(FloatMatrix data, List<String> Headers, String CSV_FILE_NAME)  {
        File csvOutputFile = new File(CSV_FILE_NAME);
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
            pw.println(String.join(",", Headers));
            pw.println(data.toString("%f", "", "",", ", "\n"));
        } catch(IOException error) {
            System.out.println("Could not save CSV.");
            error.printStackTrace();
        }
    }

    public static int getBins(FloatMatrix A, float width){
        return (int) ((A.max()-A.min())/width);
    }

    public static String getOrdinal(int i) {
        String[] sufixes = new String[] { "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th" };
        switch (i % 100) {
            case 11:
            case 12:
            case 13:
                return i + "th";
            default:
                return i + sufixes[i % 10];
        }
    }

    public static String getTitleHist(int i){
        return(getOrdinal(i) + "-" + getOrdinal(i + 1));
    }

    public static String getTitlePlot(int orders){
        StringBuilder title = new StringBuilder();
        for(int i = 0; i < orders; i++){
            title.append(getOrdinal(i)).append("\n");
        }
        return title.toString();
    }

    public static FloatMatrix cleanup(final FloatMatrix A, int neighbours, float distance, int coreCount){
        List<Integer> test = Collections.synchronizedList(new ArrayList<>());

        final FloatMatrix X = A.getColumn(3);
        final FloatMatrix Y = A.getColumn(4);

        AtomicInteger ai = new AtomicInteger(0);
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

        int[] indices =  test.stream().mapToInt(i->i).toArray();
        return A.getRows(indices);
    }

    public static FloatMatrix connectOrders(FloatMatrix intermediate, int orders, int orderCollumns){
        boolean checkMoreOrders = true;
        List<Integer> toKeep = new ArrayList<>();
        for(int row = 0; row < intermediate.rows; row ++){
            float index = intermediate.get(row, 2);

            if(intermediate.getColumn(6).eq(index).sum() == 0.0f) { //Start of chain
                toKeep.add(row);

                int connected_to = (int) intermediate.get(row, 6);
                int[] connected_indices = intermediate.getColumn(2).eq(connected_to).findIndices();

                if(connected_indices.length > 0){
                    //We already did the 0-1 connection

                    intermediate.putRow(row, recursiveSearch(intermediate, intermediate.getRow(row), connected_indices, 2, orders, orderCollumns));

                    if(checkMoreOrders){
                        connected_to = (int) intermediate.get(row, (orders-1) * orderCollumns);

                        if(intermediate.getColumn(2).eq(connected_to).sum() > 1.0f) {
                            System.out.println("There seem to be more orders than " + orders);
                            checkMoreOrders = false;
                        }
                    }

                }
            }
        }

        int[] indices = toKeep.stream().mapToInt(i->i).toArray();
        return intermediate.getRows(indices);
    }

    private static FloatMatrix recursiveSearch(FloatMatrix intermediate, FloatMatrix rowData, int[] connected_indices, int order, int max_order, int orderCollumns) {
        //i.e. 12: id, 13: x, 14: y, 15: intensity, 16: distance, 17: angle
        int[] target_range = new int[]{order * orderCollumns, 1 + (order * orderCollumns), 2 + (order * orderCollumns), 3 + (order * orderCollumns), 4 + (order * orderCollumns), 5 + (order * orderCollumns)};

        if (connected_indices.length == 1) {
            FloatMatrix target_row = intermediate.getRow(connected_indices[0]);
            rowData.put(target_range, target_row.getColumns(new int[]{6, 7, 8, 9, 10, 11}));
        } else {
            FloatMatrix intermediateResults = new FloatMatrix(connected_indices.length, rowData.columns);

            for (int i = 0; i < connected_indices.length; i++) {

                FloatMatrix target_row = intermediate.getRow(connected_indices[i]);
                intermediateResults.putRow(i, rowData.dup().put(target_range, target_row.getColumns(new int[]{6, 7, 8, 9, 10, 11})));
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

    public static int[] FilterbyLength(FloatMatrix A, int orders){
        for(int i = orders - 1; i > 1; i--){
            FloatMatrix IdCollumn = A.getColumn(i * 5);
            if(IdCollumn.sum() > 0.0f){
                return IdCollumn.ne(0.0f).findIndices();
            }
        }
        return new int[] {0};
    }

    public static int getClosestIndex(FloatMatrix A, float v){
        FloatMatrix minimised = abs(A.sub(v));
        return minimised.eq(minimised.min()).findIndices()[0];
    }

    public static FloatMatrix abs(FloatMatrix A){
        for(int i = 0; i < A.length; i++) A.put(i, Math.abs(A.get(i)));
        return A;
    }

    public static ImagePlus getImageFromPoints(FloatMatrix A, float[] reduction, int width, int height){
        float[] data = new float[width * height];
        float max_intensity = A.getColumn(2).max();
        FloatMatrix intensities = A.getColumn(2).div(max_intensity);
        for(int i = 0; i < A.rows; i++){
            data[(int)(A.get(i, 0)/reduction[0]) + (int)(A.get(i, 1)/reduction[1]) * height] = intensities.get(i);
        }

        return new ImagePlus("Points", new FloatProcessor(width, height, data, null));
    }

    public static < T extends RealType< T > > double getSum( final Iterable< T > iterable ) {
        final RealSum sum = new RealSum();

        for ( final T type : iterable )
            sum.add( type.getRealDouble() );

        return sum.getSum();
    }

    public static void addLutLegend(Plot plot, OwnColorTable ct, String label, int width, double start, double end){
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
        long sum = 0L;

        for (long l : rx) {
            sum += l;
        }

        return sum;
    }

    public static int sum(boolean[] rx){
        int sum = 0;
        for(boolean b : rx){
            sum += b ? 1 : 0;
        }
        return sum;
    }

    public static int getThresholdBin(float threshold, long[] hist){
        int bins = hist.length;
        long elementSum = sum(hist);

        threshold *= elementSum;

        long temp = 0L;
        for(int i = 0; i < bins; i++){
            temp += hist[i];
            if(temp > threshold) return i;
        }

        return 0;
    }

    public static float[][] sortMultiple(float[] a, float[] b){
        float[] tempA = a.clone();
        float[] tempB = new float[b.length];
        //float[] tempC = new float[c.length];

        Arrays.sort(tempA);
        reverse(tempA);

        for(int i = 0; i < a.length; i++){
            for(int j = 0; j < a.length; j++){
                if(a[j]==tempA[i]){
                    tempB[i] = b[j];
                    //tempC[i] = c[j];
                    break;
                }
            }
        }
        //return new float[][]{tempA, tempB, tempC};
        return new float[][]{tempA, tempB};
    }

    public static void reverse(float[] input) {
        int last = input.length - 1;
        int middle = input.length / 2;
        for (int i = 0; i <= middle; i++) {
            float temp = input[i];
            input[i] = input[last - i];
            input[last - i] = temp;
        }
    }

    public static boolean[][] checkForRetry(FloatMatrix A){
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
}
