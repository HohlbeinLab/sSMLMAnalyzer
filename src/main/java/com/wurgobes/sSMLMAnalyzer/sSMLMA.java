package com.wurgobes.sSMLMAnalyzer;


import ij.*;
import ij.plugin.PlugIn;
import org.jblas.exceptions.LapackException;
import org.scijava.plugin.Parameter;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jblas.FloatMatrix;

import static com.wurgobes.sSMLMAnalyzer.Util.*;
import static com.wurgobes.sSMLMAnalyzer.levenshtein.getTheClosestMatch;



public class sSMLMA implements PlugIn {

    private static final String COMMA_DELIMITER = ",";


    private final int[] unit_decades = {0, 0, 0, -2, -3, -6, -9, -10, -12, -15};

    private final String[] possible_options = {"id", "frame", "x", "y", "z", "intensity", "offset", "bkgstd", "sigma1", "sigma2", "uncertainty", "detections", "chi"};
    private final String[] unit_prefixes = {"null", "photons", "m", "cm", "mm", "um", "nm", "ang", "pm", "fm"};

    @Parameter
    private final String filePath = "F:\\ThesisData\\output\\output1_merged_filter2.csv";

    @Parameter
    private static final String CSV_FILE_NAME = "C:\\Users\\Martijn\\Desktop\\Thesis2020\\SpectralData\\own_output.csv";

    @Parameter
    private final boolean saveSCV = true;

    @Parameter
    private final float[] angRange = {(float) (-0.04 * Math.PI), (float)(0.04 * Math.PI)};

    @Parameter
    private final float[] distRange = {1700, 2500}; //default was {2500, 4500}


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

        int id = 0;

        FloatMatrix finalPossibilities = new FloatMatrix(0, 10);
        for(int frame = 1; frame <= frames; frame++){
            IJ.showProgress(frame, frames);

            int[] frameIndicices = data.getColumn(0).eq(frame).findIndices();

            FloatMatrix frameData = data.getRows(frameIndicices);

            FloatMatrix subtractedX = makeSubstractedMatrix(frameData.getColumn(1));
            FloatMatrix subtractedY = makeSubstractedMatrix(frameData.getColumn(2));

            FloatMatrix distances = Distance(subtractedX, subtractedY);
            FloatMatrix angles = atan2(subtractedX, subtractedY);


            int[] correctAngleAndDistance = distances.gt(distRange[0]).and(distances.lt(distRange[1]))
                    .and(angles.gt(angRange[0])).and(angles.lt(angRange[1])).findIndices();


            if(correctAngleAndDistance.length > 1){

                FloatMatrix possibilities = new FloatMatrix(correctAngleAndDistance.length, 10);

                for(int i = 0; i < correctAngleAndDistance.length; i++) {
                    int index = correctAngleAndDistance[i];
                    possibilities.putRow(i, new FloatMatrix(1, 10,
                            id++,                                                           //0
                            frame,                                                          //1
                            Math.floorDiv(index, distances.rows),                           //2
                            frameData.get(index / distances.rows, 1),    //3
                            frameData.get(index / distances.rows, 2),    //4
                            index % distances.rows,                                         //5
                            frameData.get(index % distances.rows, 1),    //6
                            frameData.get(index % distances.rows, 2),    //7
                            distances.get(index),                                           //8
                            angles.get(index)                                               //9
                    ));
                }
                FloatMatrix connectedIds = possibilities.getColumn(5);
                List<Integer> toSkip = new ArrayList<>();
                for(int i = 0; i < possibilities.rows; i++){
                    if(!toSkip.contains(i)){
                        FloatMatrix identicalIds = connectedIds.eq(connectedIds.get(i));
                        float identicalIdSum = identicalIds.sum();
                        if(identicalIdSum > 1.0f){

                            float cumX = possibilities.get(i, 3);
                            float cumY = possibilities.get(i, 4);

                            for(int index2 : identicalIds.findIndices()) {
                                if(index2 != i) {
                                    cumX += possibilities.get(index2, 3);
                                    cumY += possibilities.get(index2, 4);
                                    toSkip.add(index2);
                                }
                            }

                            FloatMatrix temp_arr = possibilities.getRow(i).dup();
                            temp_arr.put(3, cumX/identicalIdSum);
                            temp_arr.put(4, cumY/identicalIdSum);

                            finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, temp_arr);

                        } else {
                            finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, possibilities.getRow(i));
                        }
                    }
                }

            } else if(correctAngleAndDistance.length > 0){
                int index = correctAngleAndDistance[0];
                finalPossibilities = FloatMatrix.concatVertically(finalPossibilities, new FloatMatrix(1, 10,
                        id++,
                        frame,
                        Math.floorDiv(index, distances.rows),
                        frameData.get( index/distances.rows, 1),
                        frameData.get( index/distances.rows, 2),
                        index%distances.rows,
                        frameData.get( index%distances.rows, 1),
                        frameData.get( index%distances.rows, 2),
                        distances.get(index),
                        angles.get(index)));
            }
        }

        //Fix ids
        id = 0;
        String[] Header = {
                "id","frame",
                "start index",
                "start x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]",
                "start y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]",
                "end index",
                "end x [" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "]",
                "end y [" + unit_prefixes[unitsIndices[revOptionsIndices[3]]] + "]",
                "distance [" + (unit_prefixes[unitsIndices[revOptionsIndices[2]]].equals(unit_prefixes[unitsIndices[revOptionsIndices[3]]]) ? unit_prefixes[unitsIndices[revOptionsIndices[2]]]:("(" + unit_prefixes[unitsIndices[revOptionsIndices[2]]] + "*" +unit_prefixes[unitsIndices[revOptionsIndices[2]]] + ")^½"))+ "]",
                "angle"
        };
        System.out.println(Arrays.toString(Header));
        for(int i = 0; i < finalPossibilities.rows; i++){
            finalPossibilities.put(i, 0, id++);
            System.out.println(finalPossibilities.getRow(i));
        }

        processingTime = System.nanoTime() - processingTime;
        System.out.println("\nProcessing data took " + String.format("%.3f", processingTime/1000000000) + " s");



        if(saveSCV) SaveCSV(finalPossibilities, Header);

    }

    public void SaveCSV(FloatMatrix data, String[] Headers)  {
        File csvOutputFile = new File(CSV_FILE_NAME);
        try (PrintWriter pw = new PrintWriter(csvOutputFile)) {
                pw.println(String.join(",", Headers));
                pw.println(data.toString("%f", "", "",", ", "\n"));
        } catch(IOException error) {
            System.out.println("Could not save CSV.");
            error.printStackTrace();
        }

    }



    public static void main(String[] args) {
        new ImageJ();
        IJ.runPlugIn(sSMLMA.class.getName(), "");
        System.exit(0);
    }
}

