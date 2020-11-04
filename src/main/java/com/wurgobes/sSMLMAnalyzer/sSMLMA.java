package com.wurgobes.sSMLMAnalyzer;


import ij.*;
import ij.plugin.PlugIn;
import org.jblas.exceptions.LapackException;
import org.scijava.plugin.Parameter;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jblas.FloatMatrix;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;
import static org.jblas.Geometry.pairwiseSquaredDistances;


public class sSMLMA implements PlugIn {

    private static final String COMMA_DELIMITER = ",";

    private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};

    @Parameter
    private String filePath = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\small.csv";

    @Parameter
    private final boolean saveSCV = false;

    @Parameter
    private final double[] angRange = {-0.04 * Math.PI, 0.04 * Math.PI};

    @Parameter
    private final int[] distRange = {2500, 4500};


    private FloatMatrix getRecordFromLine(String line, int size) {
        FloatMatrix values = new FloatMatrix(1, size);
        try (Scanner rowScanner = new Scanner(line)) {
            rowScanner.useDelimiter(COMMA_DELIMITER);
            int i = 0;
            while (rowScanner.hasNextFloat()) {
                values.put(i, rowScanner.nextFloat());
                i++;
            }
        }
        return values;
    }

    @Override
    public void run(String arg) {
        double csvTime = System.nanoTime();

        FloatMatrix floatMatrix = null;
        List<String> collumns = new ArrayList<>();

        OwnFloatMatrix ownFloatMatrix = new OwnFloatMatrix();

        try {
            floatMatrix =  ownFloatMatrix.loadCSVFile(filePath);
            collumns = ownFloatMatrix.collumns;
        } catch (IOException e) {
            System.out.println("File not found.");
        } catch (LapackException e) {
            e.printStackTrace();
        } catch (Exception e){
            System.out.println("Does the csv start with a header?");
        }

        assert !collumns.isEmpty();
        assert floatMatrix != null;


        int[] revOptionsIndices = new int[possible_options.length];
        int[] unitsIndices = new int[collumns.size()];


        Pattern pattern = Pattern.compile("(\\w+)( [ (\\[](\\w+)[)\\] ])?");

        for(int i = 0; i < collumns.size(); i ++){
            String header = collumns.get(i);
            Matcher matcher = pattern.matcher(header);
            if(matcher.find()){
                revOptionsIndices[getTheClosestMatch(possible_options, matcher.group(1))] = i;
                unitsIndices[i] = getTheClosestMatch(unit_prefixes, matcher.group(3));
            }
        }

        final FloatMatrix data = floatMatrix.getColumns(new int[]{revOptionsIndices[1], revOptionsIndices[2], revOptionsIndices[3]});

        final int frames = (int) data.getColumn(0).max(); //1 indexed

        csvTime = System.nanoTime() - csvTime;
        System.out.println("Loading CSV took " + String.format("%.3f", csvTime/1000000000) + " s");
        double processingTime = System.nanoTime();


        for(int frame = 1; frame <= frames; frame++){
            IJ.showProgress(frame, frames);

            int[] frameIndicices = data.getColumn(0).eq(frame).findIndices();

            FloatMatrix frameData = data.getRows(frameIndicices);

            FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
            FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

            FloatMatrix distances = Distances(frameData.getColumn(1), frameData.getColumn(2));
            FloatMatrix angles = AnglesBetweenPoints(frameData.getColumn(1), frameData.getColumn(2));

            System.out.println(angles.getColumn(0));
            System.exit(0);
        }

        processingTime = System.nanoTime() - processingTime;
        System.out.println("Processing data took " + String.format("%.3f", processingTime/1000000000) + " s");




    }

    public static void main(String[] args) {
        new ImageJ();
        IJ.runPlugIn(sSMLMA.class.getName(), "");
        System.exit(0);
    }
}

