package com.wurgobes.sSMLMAnalyzer;

import org.jblas.FloatMatrix;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

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
        FloatMatrix domainTwo = atan(intermediate.sub(x).div(y)).mul(2);

        for(int index : x.gt(0).findIndices()) result.put(index, domainOne.get(index));
        for(int index : x.le(0).and(y.ne(0)).findIndices()) result.put(index, domainTwo.get(index));
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

    public static void mirror(FloatMatrix A, int col){
        FloatMatrix temp = A.getColumn(col);
        float mean = temp.mean();
        temp.muli(-1.0f).addi(mean * 2);
        A.putColumn(col, temp);
    }
}
