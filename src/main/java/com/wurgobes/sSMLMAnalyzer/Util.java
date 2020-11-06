package com.wurgobes.sSMLMAnalyzer;

import org.jblas.FloatMatrix;


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
        for(int index : x.eq(0).and(y.eq(0)).findIndices()) result.put(index, Float.NaN);

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
            subtracted.putRow(point, A.sub(a));
        }
        return subtracted;
    }

    public static double[] toDouble(FloatMatrix A){
        double[] result = new double[A.length];

        for(int i = 0; i < A.length; i++) result[i] = (double) A.get(i);

        return result;
    }
}
