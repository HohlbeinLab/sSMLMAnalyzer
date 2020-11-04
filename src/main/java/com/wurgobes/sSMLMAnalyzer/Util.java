package com.wurgobes.sSMLMAnalyzer;

import org.jblas.FloatMatrix;
import org.jblas.*;

import java.util.Arrays;

import static org.jblas.MatrixFunctions.sqrt;
import static org.jblas.MatrixFunctions.atan;

//Very much not a perfect solution
//https://en.wikipedia.org/wiki/Atan2

public class Util {

    public static FloatMatrix atan2(FloatMatrix x, FloatMatrix y){
        FloatMatrix result = new FloatMatrix(x.getRows());

        FloatMatrix intermediate = sqrt(x.mul(x).add(y.mul(y)));

        FloatMatrix domainOne = atan(y.div(intermediate.add(x))).mul(2);
        FloatMatrix domainTwo = atan(intermediate.sub(x).div(y)).mul(2);

        for(int index : x.gt(0).findIndices()) result.put(index, domainOne.get(index));
        for(int index : x.le(0).and(y.ne(0)).findIndices()) result.put(index, domainTwo.get(index));
        for(int index : x.lt(0).and(y.eq(0)).findIndices()) result.put(index, (float) Math.PI);
        for(int index : x.eq(0).and(y.eq(0)).findIndices()) result.put(index, Float.NaN);

        return result;
    }

    public static FloatMatrix AnglesBetweenPoints(FloatMatrix X, FloatMatrix Y){
        if (X.rows != Y.rows)
            throw new IllegalArgumentException(
                    "Matrices must have same number of rows");

        int dataPoints = X.rows;
        FloatMatrix angles = new FloatMatrix(dataPoints, dataPoints);

        FloatMatrix anglesX = new FloatMatrix(dataPoints, dataPoints);
        FloatMatrix anglesY = new FloatMatrix(dataPoints, dataPoints);


        //System.out.println(angles);
        //System.out.println(angles.getRows());
        //System.out.println(angles.getColumns());

        for(int point = 0; point < dataPoints; point++){
            float x = X.get(point);
            float y = Y.get(point);

            FloatMatrix temp_x = X.sub(x);
            FloatMatrix temp_y = Y.sub(y);

            angles.putColumn(point, atan2(temp_x, temp_y));

            anglesX.putColumn(point, X.sub(x));
            anglesY.putColumn(point, Y.sub(y));
        }


        return angles;
    }

    public static FloatMatrix Distance(FloatMatrix X, FloatMatrix Y){
        return sqrt(X.mul(X).addi(Y.mul(Y)));
    }

    public static FloatMatrix Distances(FloatMatrix X, FloatMatrix Y) {
        if (X.rows != Y.rows)
            throw new IllegalArgumentException(
                    "Matrices must have same number of rows");

        int dataPoints = X.rows;

        FloatMatrix distancesX = new FloatMatrix(dataPoints, dataPoints);
        FloatMatrix distancesY = new FloatMatrix(dataPoints, dataPoints);

        for(int point = 0; point < dataPoints; point++){
            float x = X.get(point);
            float y = Y.get(point);

            distancesX.putColumn(point, X.sub(x));
            distancesY.putColumn(point, Y.sub(y));
        }



        return Distance(distancesX, distancesY);
    }

    public static FloatMatrix makeSubstractedMatrix(FloatMatrix A){
        if(A.columns > 1)
            throw new IllegalArgumentException("Matrix can have only 1 collumn");

        int dataPoints = A.rows;

        FloatMatrix subtracted = new FloatMatrix(dataPoints, dataPoints);

        for(int point = 0; point < dataPoints; point++){
            float a = A.get(point);
            subtracted.putColumn(point, A.sub(a));

        }
        return subtracted;
    }
}
