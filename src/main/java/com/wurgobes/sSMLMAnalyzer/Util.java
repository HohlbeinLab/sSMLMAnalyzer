package com.wurgobes.sSMLMAnalyzer;

import ij.IJ;
import org.jblas.FloatMatrix;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static ij.plugin.filter.MaximumFinder.findMaxima;
import static org.jblas.MatrixFunctions.sqrt;
import static org.jblas.MatrixFunctions.atan;



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

    public static int[] getMaxima(int[] values, int maxRatio){
        float maxValue = (float) Arrays.stream(values).max().getAsInt();
        return findMaxima(toDouble(values), maxValue/maxRatio, 0);
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
}
