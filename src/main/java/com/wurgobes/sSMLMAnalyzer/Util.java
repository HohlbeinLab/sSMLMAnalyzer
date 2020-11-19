package com.wurgobes.sSMLMAnalyzer;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealPointSampleList;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.type.numeric.real.FloatType;
import org.jblas.FloatMatrix;


import javax.xml.stream.FactoryConfigurationError;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static ij.plugin.filter.MaximumFinder.findMaxima;
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

        for(int i = 0; i < A.length; i++) result[i] = (double) A.get(i);

        return result;
    }

    public static double[] toDouble(int[] A){
        double[] result = new double[A.length];

        for(int i = 0; i < A.length; i++) result[i] = A[i];

        return result;
    }

    public static double[] toDouble(float v){
        return new double[] {v};
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
            title.append(getOrdinal(i) + "\n");
        }
        return title.toString();
    }

    public static FloatMatrix cleanup(FloatMatrix A, int neighbours, float distance){
        List<Integer> indices = new ArrayList<>();

        FloatMatrix X = A.getColumn(3);
        FloatMatrix Y = A.getColumn(4);


        for(int i = 0; i < A.rows; i++){
            IJ.showProgress(i, A.rows);
            float x = A.get(i, 3);
            float y = A.get(i, 4);

            FloatMatrix distances = Distance(X.sub(x), Y.sub(y));

            if(distances.lt(distance).sum() > neighbours + 1) indices.add(i);
        }

        int[] indices2 =  indices.stream().mapToInt(i->i).toArray();
        return A.getRows(indices2);
    }

    public static FloatMatrix connectOrders(FloatMatrix intermediate, int orders, int orderCollumns){
        boolean checkMoreOrders = true;
        List<Integer> toKeep = new ArrayList<>();
        for(int row = 0; row < intermediate.rows; row ++){
            float index = intermediate.get(row, 2);

            if(intermediate.getColumn(5).eq(index).sum() == 0.0f) { //Start of chain
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
        //i.e. 10: id, 11: x, 12: y, 13: distance, 14: angle
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

    public static RealPointSampleList<FloatType> createPointList(FloatMatrix floatMatrix, int[] xyz, int intensity, float reduce, float minZ) {
        // the number of dimensions
        final int numDimensions = xyz.length;


        // a list of Samples with coordinates
        RealPointSampleList< FloatType > elements = new RealPointSampleList<>( numDimensions );

        for ( int i = 0; i < floatMatrix.rows; ++i ) {
            RealPoint point = new RealPoint( numDimensions );

            for ( int d = 0; d < numDimensions - 1; ++d )
                point.setPosition( floatMatrix.get(i, xyz[d])/reduce, d );

            if(numDimensions > 2)
                point.setPosition(floatMatrix.get(i, xyz[numDimensions - 1])-minZ, numDimensions - 1);

            float value = intensity == -1 ? 1f : floatMatrix.get(i, intensity);
            elements.add(point, new FloatType(value));
        }

        return elements;
    }

    public static long[] toLong(float[] arr){
        long[] result = new long[arr.length];
        for(int i = 0;i < result.length; i ++)
            result[i] = (long) arr[i];

        return result;
    }

    public static FloatProcessor getImageFromPoints(FloatMatrix A, int[] reduction, int width, int height){
        float[] data = new float[width * height];
        float max_intensity = A.getColumn(2).max();
        FloatMatrix intensities = A.getColumn(2).div(max_intensity);
        for(int i = 0; i < A.rows; i++){
            data[(int)(A.get(i, 0)/reduction[0]) + (int)(A.get(i, 1)/reduction[1]) * height] = intensities.get(i);
        }

        return new FloatProcessor(width, height, data, null);
    }
}
