package com.wurgobes.sSMLMAnalyzer;

/*
Spectral Super Resolution Pair Finder
(c) 2021 Martijn Gobes, Wageningen University.

This file contains a loader that takes a csv file with a header
The Header is parsed into an arraylist and the values into a jblas Floatmatrix
Default and fallback delimiter is a comma

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

import org.jblas.FloatMatrix;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Pattern;


public class OwnFloatMatrixLoader {


    private List<String> columns = new ArrayList<>();

    public FloatMatrix loadCSVFile(String filename) throws IOException {


        BufferedReader is = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));

        List<FloatMatrix> rows = new LinkedList<>();
        String line = is.readLine();
        int columns = -1;

        String delimiter = ",";
        if(line.contains(";")){
            delimiter = ";";
        } else if(line.contains("\t")) {
            delimiter = "\t";
        }

        final Pattern DELIMITER_PATTERN  = Pattern.compile(";");
        final String DELIMITER = delimiter;

        this.columns = Arrays.asList(line.split(DELIMITER));


        while ((line = is.readLine()) != null) {
            String[] elements = DELIMITER_PATTERN.split(line);
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
                row.put(c, Float.parseFloat(elements[c]));
            }
            rows.add(row);
        }
        is.close();

        System.out.println("Done reading file: " + filename);

        FloatMatrix result = new FloatMatrix(rows.size(), columns);
        int r = 0;
        for (FloatMatrix row : rows) {
            result.putRow(r, row);
            r++;
        }
        return result;
    }

    public List<String> getColumns() {return columns;}
}
