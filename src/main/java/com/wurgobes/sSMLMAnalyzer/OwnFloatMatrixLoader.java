package com.wurgobes.sSMLMAnalyzer;

import org.jblas.FloatMatrix;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Pattern;


public class OwnFloatMatrixLoader {
    static final Pattern COMMA = Pattern.compile(",");
    static final String COMMA_DELIMITER = ",";

    public List<String> collumns = new ArrayList<>();

    public FloatMatrix loadCSVFile(String filename) throws IOException {
        BufferedReader is = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));

        List<FloatMatrix> rows = new LinkedList<>();
        String line;
        int columns = -1;


        collumns = Arrays.asList(is.readLine().split(COMMA_DELIMITER));


        while ((line = is.readLine()) != null) {
            String[] elements = COMMA.split(line);
            int numElements = elements.length;
            if (elements[0].length() == 0) {
                numElements--;
            }
            if (elements[elements.length - 1].length() == 0) {
                numElements--;
            }

            if (columns == -1) {
                columns = numElements;
            } else {
                if (columns != numElements) {
                    throw new IOException("Number of elements changes in line " + line + ".");
                }
            }

            FloatMatrix row = new FloatMatrix(columns);
            for (int c = 0; c < columns; c++) {
                row.put(c, Float.valueOf(elements[c]));
            }
            rows.add(row);
        }
        is.close();

        System.out.println("Done reading file");

        FloatMatrix result = new FloatMatrix(rows.size(), columns);
        int r = 0;
        for (FloatMatrix row : rows) {
            result.putRow(r, row);
            r++;
        }
        return result;
    }
}
